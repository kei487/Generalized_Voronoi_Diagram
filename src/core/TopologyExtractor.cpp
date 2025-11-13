#include "gvd_topo/core/TopologyExtractor.hpp"
#include "gvd_topo/core/OccupancyGrid.hpp"
#include <cstdint>
#include <queue>
#include <cmath>
#include <limits>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <vector>
#include <unordered_set>

#include <iostream>

#ifdef GVD_TOPO_WITH_OPENCV
#include <opencv2/imgproc.hpp>
#include <opencv2/core.hpp>
#endif

namespace gvd_topo {

TopologyExtractor::TopologyExtractor() = default;
TopologyExtractor::TopologyExtractor(const Params& p) : params_(p) {}

static inline int idx(int x, int y, int w) { return y * w + x; }

// 8-neighborhood offsets
static const int dx8[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
static const int dy8[8] = { -1,-1,-1,  0, 0,  1, 1, 1 };

// Node pixel structure
struct NodePix { int x; int y; };

// Temporary node structure
struct TempNode {
    int id;
    int x, y;  // pixel coordinates
    double world_x, world_y;
};

// ============================================================================
// Helper Functions for Topology Extraction
// ============================================================================

// Optimized Zhang-Suen thinning algorithm with parallel processing
static std::vector<uint8_t> zhangSuenThinning(const std::vector<uint8_t>& input, int width, int height) {
    std::vector<uint8_t> img = input;
    std::vector<uint8_t> marker(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
    
    bool changed = true;
    int iteration = 0;
    // Adjust max iterations based on map size
    const int max_iterations = 100;//(width * height > 10000000) ? 20 : 100;
    
    while (changed && iteration < max_iterations) {
        changed = false;
        
        // Sub-iteration 1: Mark pixels to delete
    #ifdef GVD_TOPO_WITH_OPENMP
        #pragma omp parallel for schedule(dynamic, 64)
    #endif
        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                int idx_center = idx(x, y, width);
                if (img[idx_center] == 0) continue;
                
                // Get 8 neighbors (inline for speed)
                bool P2 = img[idx(x, y-1, width)] > 0;  // top
                bool P3 = img[idx(x+1, y-1, width)] > 0;  // top-right
                bool P4 = img[idx(x+1, y, width)] > 0;  // right
                bool P5 = img[idx(x+1, y+1, width)] > 0;  // bottom-right
                bool P6 = img[idx(x, y+1, width)] > 0;  // bottom
                bool P7 = img[idx(x-1, y+1, width)] > 0;  // bottom-left
                bool P8 = img[idx(x-1, y, width)] > 0;  // left
                bool P9 = img[idx(x-1, y-1, width)] > 0;  // top-left
                
                // Count non-zero neighbors
                int B = P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9;
                
                // Count 0->1 transitions
                int A = (!P2 && P3) + (!P3 && P4) + (!P4 && P5) + (!P5 && P6) +
                        (!P6 && P7) + (!P7 && P8) + (!P8 && P9) + (!P9 && P2);
                
                // Conditions for sub-iteration 1
                if (B >= 2 && B <= 6 && A == 1 && !(P2 && P4 && P6) && !(P4 && P6 && P8)) {
                    marker[idx_center] = 1;
                }
            }
        }
        
        // Delete marked pixels
        const size_t total_pixels = static_cast<size_t>(width) * static_cast<size_t>(height);
        for (size_t i = 0; i < total_pixels; ++i) {
            if (marker[i] == 1) {
                img[i] = 0;
                marker[i] = 0;
                changed = true;
            }
        }
        
        // Sub-iteration 2: Mark pixels to delete
        #ifdef GVD_TOPO_WITH_OPENMP
        #pragma omp parallel for schedule(dynamic, 64)
        #endif
        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                int idx_center = idx(x, y, width);
                if (img[idx_center] == 0) continue;
                
                // Get 8 neighbors (inline for speed)
                bool P2 = img[idx(x, y-1, width)] > 0;
                bool P3 = img[idx(x+1, y-1, width)] > 0;
                bool P4 = img[idx(x+1, y, width)] > 0;
                bool P5 = img[idx(x+1, y+1, width)] > 0;
                bool P6 = img[idx(x, y+1, width)] > 0;
                bool P7 = img[idx(x-1, y+1, width)] > 0;
                bool P8 = img[idx(x-1, y, width)] > 0;
                bool P9 = img[idx(x-1, y-1, width)] > 0;
                
                // Count non-zero neighbors
                int B = P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9;
                
                // Count 0->1 transitions
                int A = (!P2 && P3) + (!P3 && P4) + (!P4 && P5) + (!P5 && P6) +
                        (!P6 && P7) + (!P7 && P8) + (!P8 && P9) + (!P9 && P2);
                
                // Conditions for sub-iteration 2
                if (B >= 2 && B <= 6 && A == 1 && !(P2 && P4 && P8) && !(P2 && P6 && P8)) {
                    marker[idx_center] = 1;
                }
            }
        }
        
        // Delete marked pixels
        for (size_t i = 0; i < total_pixels; ++i) {
            if (marker[i] == 1) {
                img[i] = 0;
                marker[i] = 0;
                changed = true;
            }
        }
        
        ++iteration;
    }
    
    return img;
}


// Extract nodes from Hough transform line detection (alternative to degree-based method)
static TopologicalMap extractNodesFromHoughLines(const std::vector<uint8_t>& is_skel, int width, int height,
     double resolution,double min_length) {
#ifdef GVD_TOPO_WITH_OPENCV
    std::vector<TopoNode> nodes;
    std::vector<TopoEdge> edges;

    // Convert skeleton to OpenCV Mat
    cv::Mat skel_img(height, width, CV_8UC1);
    for (int y = 0; y < height; ++y) {
        uint8_t* row = skel_img.ptr<uint8_t>(y);
        for (int x = 0; x < width; ++x) {
            row[x] = is_skel[idx(x, y, width)];
        }
    }
    
    // Apply Hough Line Transform (Probabilistic)
    std::vector<cv::Vec4i> lines;
    const double rho = 1.0;              // Distance resolution in pixels
    const double theta = CV_PI / 180.0;  // Angular resolution in radians (1 degree)
    const int threshold = 5;             // Minimum number of votes (reduced for sensitivity)
    const double min_line_length = 5.0;  // Minimum line length (reduced)
    const double max_line_gap = 10.0;    // Maximum gap between line segments (increased)
    
    cv::HoughLinesP(skel_img, lines, rho, theta, threshold, min_line_length, max_line_gap);
    
    // Debug: print number of detected lines (temporary)
    // std::cout << "Hough transform detected " << lines.size() << " lines" << std::endl;
    
    // Extract endpoints and intersections from detected lines
    std::set<std::pair<int, int>> node_set; // Use set to avoid duplicates
    TopoNode tmp_node;
    TopoEdge tmp_edge;
    int tmp_id = 0, node_id = 0;   

    // Add endpoints
    for (const auto& line : lines) {
        int x1 = line[0], y1 = line[1];
        int x2 = line[2], y2 = line[3];
        
        // Calculate length in pixels, then convert to meters
        double length_px = std::hypot(x2 - x1, y2 - y1);
        tmp_edge.length = length_px * resolution;

        if(length_px < min_length) {
            continue;
        }

        tmp_node.id = node_id++;
        tmp_node.x = x1;
        tmp_node.y = y1;
        //tmp_node.edge_ids.push_back(tmp_id);
        nodes.push_back(tmp_node);
        
        tmp_node.id = node_id++;
        tmp_node.x = x2;
        tmp_node.y = y2;
        //tmp_node.edge_ids.push_back(tmp_id);
        nodes.push_back(tmp_node);
    
        tmp_edge.id = tmp_id++;
        tmp_edge.u = node_id - 2;
        tmp_edge.v = node_id - 1;
        
        edges.push_back(tmp_edge);


    }

    TopologicalMap topo;
    topo.nodes = nodes;
    topo.edges = edges;
    return topo;
#else
    (void)is_skel; (void)width; (void)height;
    return TopologicalMap(); // Return empty if OpenCV not available
        #endif
}

// Bresenham's line algorithm to get all pixels on a line
static std::vector<std::pair<int, int>> bresenhamLine(int x0, int y0, int x1, int y1) {
    std::vector<std::pair<int, int>> pixels;
    
    int dx = std::abs(x1 - x0);
    int dy = std::abs(y1 - y0);
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx - dy;
    
    int x = x0;
    int y = y0;
    
    while (true) {
        pixels.push_back({x, y});
        
        if (x == x1 && y == y1) break;
        
        int e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x += sx;
        }
        if (e2 < dx) {
            err += dx;
            y += sy;
        }
    }
    
    return pixels;
}

// Check if a line path is collision-free on the occupancy grid
static bool isPathClearOnGrid(const OccupancyGrid& grid, int x0, int y0, int x1, int y1) {
    auto pixels = bresenhamLine(x0, y0, x1, y1);
    
    for (const auto& [x, y] : pixels) {
        // Check bounds
        if (!grid.inBounds(x, y)) {
            return false;
        }
        
        // Check if cell is occupied
        int8_t cell_value = grid.data[grid.index(x, y)];
        // Consider occupied (100) and unknown (-1) as obstacles
        if (cell_value != static_cast<int8_t>(Cell::Free)) {
            return false;
        }
    }
    
    return true;
}

// Connect isolated endpoint nodes by finding closest collision-free partners
static TopologicalMap connectEndpoints(
    const OccupancyGrid& grid, 
    const TopologicalMap& topo, 
    double max_search_radius,
    bool is_best) {
    
    TopologicalMap result = topo;

    std::map<int, int> degree_map;
    for(const auto& edge : result.edges) {
        degree_map[edge.u]++;
        degree_map[edge.v]++;
    }

    // Identify endpoints (degree == 1)
    std::vector<int> endpoint_ids;
    for (const auto& node : result.nodes) {
        if (degree_map[node.id] == 1 || degree_map[node.id] == 0) {
            endpoint_ids.push_back(node.id);
        }
    }
    std::cout << "    Debug endpoint_ids size: " << endpoint_ids.size() << std::endl;
    
    // Build spatial grid for faster neighbor search
    double max_search_radius_px = max_search_radius / grid.resolution;
    int grid_size = static_cast<int>(std::max(1.0, max_search_radius_px)) + 1;

    std::map<std::pair<int,int>, std::vector<int>> spatial_grid;
    for (auto& node : result.nodes) {
        int gx = static_cast<int>(std::round(node.x / grid.resolution));
        int gy = static_cast<int>(std::round(node.y / grid.resolution));
        spatial_grid[{gx, gy}].push_back(node.id);
    }
    
    // set next edge id
    int next_edge_id = result.edges.empty() ? 0 : 
        (*std::max_element(result.edges.begin(), result.edges.end(), 
            [](const TopoEdge& a, const TopoEdge& b) { return a.id < b.id; })).id + 1;
    
    // Set to track which endpoints have been connected
    std::set<int> connected_endpoints;
    
    //debug
    int num_edge_add=0;

    //loop through all endpoints
    for (int endpoint_id : endpoint_ids) {
        // Skip if already connected in this pass
        if (connected_endpoints.count(endpoint_id)) continue;
        
        const auto& src_node = result.nodes[endpoint_id];
        
        // Node coordinates are already in pixels (from Hough transform)
        int src_x = static_cast<int>(std::round(src_node.x));
        int src_y = static_cast<int>(std::round(src_node.y));
        int src_gx = src_x / grid.resolution;
        int src_gy = src_y / grid.resolution;
        
        
        //Set to track the best edge
        int best_candidate_id = -1;
        double best_distance_px = std::numeric_limits<double>::max();
    
       //test
        bool is_edge_add = false;

        // Check cells in the grid_size range
        for (int dgy = -grid_size/2; dgy <= grid_size/2; ++dgy) {
            for (int dgx = -grid_size/2; dgx <= grid_size/2; ++dgx) {
                auto it = spatial_grid.find({src_gx + dgx, src_gy + dgy});
                if (it == spatial_grid.end()) continue;
                for (int candidate_id : it->second) {
                    // Skip self and already connected
                    if (candidate_id == endpoint_id) {
                        // std::cout << "    Debug already connect: " 
                        // << candidate_id << ", "<< endpoint_id << std::endl;
                        continue;
                    }
                    // std::cout << "    Debug find to connect pair: " 
                    //     << candidate_id << ", "<< endpoint_id << std::endl;
                    const auto& cand_node = result.nodes[candidate_id];
                    
                    // Candidate coordinates are also in pixels
                    int cand_x = static_cast<int>(std::round(cand_node.x));
                    int cand_y = static_cast<int>(std::round(cand_node.y));
                    
                    // Check if path is collision-free
                    if (!isPathClearOnGrid(grid, src_x, src_y, cand_x, cand_y)) {
                        //std::cout << "    Debug not pass" << std::endl;
                        continue;
                    }
                    //std::cout << "    Debug pass" << std::endl;

                    //test
                    is_edge_add = true;
                    // calculate distance between source and candidate node
                    double distance_px = std::hypot(src_x - cand_x, src_y - cand_y);

                    if(is_best) {
                        // check if the edge is the best edge
                        if(distance_px < best_distance_px) {
                            best_distance_px = distance_px;
                            best_candidate_id = candidate_id;
                        }
                    }else{
                        //add edge between endpoint and candidate node
                        TopoEdge new_edge;
                        new_edge.id = next_edge_id++;
                        new_edge.u = endpoint_id;
                        new_edge.v = candidate_id;
                        // Convert pixel distance to meters
                        new_edge.length = distance_px * grid.resolution;
                    
                        result.edges.push_back(new_edge);

                        // Mark both endpoints as connected
                        connected_endpoints.insert(endpoint_id);
                        connected_endpoints.insert(candidate_id);
                    }
                    
                }
            }
        }

        if(is_best) {
            //add edge between endpoint and candidate node
            TopoEdge new_edge;
            new_edge.id = next_edge_id++;
            new_edge.u = endpoint_id;
            new_edge.v = best_candidate_id;
            // Convert pixel distance to meters
            new_edge.length = best_distance_px * grid.resolution;
        
            result.edges.push_back(new_edge);
        
            // Mark both endpoints as connected
            connected_endpoints.insert(endpoint_id);
            connected_endpoints.insert(best_candidate_id);
        }
        
        if(is_edge_add){
            num_edge_add++;
        }

    }

    std::cout << "    Debug add edgge: " << num_edge_add << std::endl;

    return result;
}

// Merge nearby nodes that are connected by short edges
static TopologicalMap mergeNearbyNodes(const TopologicalMap& topo, double merge_threshold, double resolution) {
    if (topo.nodes.empty()) return topo;
    
    // Union-Find data structure
    std::vector<int> parent(topo.nodes.size());
    for (size_t i = 0; i < parent.size(); ++i) {
        parent[i] = i;
    }
    
    auto find = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]]; // Path compression
            x = parent[x];
        }
        return x;
    };
    
    auto unite = [&](int x, int y) {
        x = find(x);
        y = find(y);
        if (x != y) {
            parent[y] = x;
        }
    };
    
    // Group nodes connected by edges shorter than threshold
    for (const auto& edge : topo.edges) {
        if (edge.length <= merge_threshold) {
            unite(edge.u, edge.v);
        }
    }
    
    // Build groups: map from root to list of node IDs
    std::map<int, std::vector<int>> groups;
    for (size_t i = 0; i < topo.nodes.size(); ++i) {
        int root = find(i);
        groups[root].push_back(i);
    }
    
    // Create merged nodes (compute centroid for each group)
    TopologicalMap result;
    std::map<int, int> old_to_new_id; // old node ID -> new node ID
    
    int new_node_id = 0;
    for (const auto& [root, group_ids] : groups) {
        TopoNode merged_node;
        merged_node.id = new_node_id++;
        
        // Compute centroid
        double sum_x = 0.0, sum_y = 0.0;
        for (int old_id : group_ids) {
            sum_x += topo.nodes[old_id].x;
            sum_y += topo.nodes[old_id].y;
        }
        merged_node.x = sum_x / group_ids.size();
        merged_node.y = sum_y / group_ids.size();
        
        result.nodes.push_back(merged_node);
        
        // Map all old IDs in this group to the new merged node ID
        for (int old_id : group_ids) {
            old_to_new_id[old_id] = merged_node.id;
        }
    }
    
    // Rebuild edges with new node IDs
    std::set<std::pair<int, int>> added_edges; // To avoid duplicate edges
    int new_edge_id = 0;
    
    for (const auto& old_edge : topo.edges) {
        int new_u = old_to_new_id[old_edge.u];
        int new_v = old_to_new_id[old_edge.v];
        
        // Skip self-loops (edges within merged group)
        if (new_u == new_v) continue;
        
        // Ensure consistent edge direction (smaller ID first)
        if (new_u > new_v) std::swap(new_u, new_v);
        
        // Skip duplicate edges
        if (added_edges.count({new_u, new_v})) continue;
        added_edges.insert({new_u, new_v});
        
        // Create new edge
        TopoEdge new_edge;
        new_edge.id = new_edge_id++;
        new_edge.u = new_u;
        new_edge.v = new_v;
        
        // Recalculate length based on merged node positions (in pixels, convert to meters)
        double dx = result.nodes[new_u].x - result.nodes[new_v].x;
        double dy = result.nodes[new_u].y - result.nodes[new_v].y;
        double length_px = std::hypot(dx, dy);
        new_edge.length = length_px * resolution;
        
        result.edges.push_back(new_edge);

        // add edge id to both nodes
        //result.nodes[new_u].edge_ids.push_back(new_edge.id);
        //result.nodes[new_v].edge_ids.push_back(new_edge.id);
    }
    
    return result;
}


// ============================================================================
// Main Topology Extraction Function
// ============================================================================

TopologicalMap TopologyExtractor::run(const OccupancyGrid& grid, const std::vector<uint8_t>& gvd_mask, int width, int height, double resolution) const {
    TopologicalMap topo;
    if (gvd_mask.empty() || width <= 0 || height <= 0) return topo;

    // Phase 0: Apply Zhang-Suen thinning as preprocessing
    std::vector<uint8_t> is_skel = zhangSuenThinning(gvd_mask, width, height);

    // Phase 1: Extract nodes and edges from Hough transform line detection
    double min_length = 5.0 * resolution;
    topo = extractNodesFromHoughLines(is_skel, width, height, resolution, min_length);

    std::cout << "   Debug topo size1: " << topo.nodes.size() << ", edges=" << topo.edges.size() << std::endl;


    // Phase 1.5: Merge nearby nodes connected by short edges
    const double merge_threshold = 10.0 * resolution;  // Merge nodes closer than 1 pixel (in meters)
    topo = mergeNearbyNodes(topo, merge_threshold, resolution);

    std::cout << "   Debug topo size1.5: " << topo.nodes.size() << ", edges=" << topo.edges.size() << std::endl;

    // Phase 2: Connect isolated endpoint nodes
    // Adjust search radius based on map size to avoid performance issues
    double max_search_radius = 10.0 * resolution;  // meters
    topo = connectEndpoints(grid, topo, max_search_radius, false);

    std::cout << "   Debug topo size2: " << topo.nodes.size() << ", edges=" << topo.edges.size() << std::endl;

    // Phase 3: Merge nearby nodes connected by short edges
    //const double merge_threshold = 10.0 * resolution;  // Merge nodes closer than 1 pixel (in meters)
    topo = mergeNearbyNodes(topo, merge_threshold, resolution);

    std::cout << "   Debug topo size3: " << topo.nodes.size() << ", edges=" << topo.edges.size() << std::endl;

    max_search_radius = 30.0 * resolution;
    topo = connectEndpoints(grid, topo, max_search_radius, false);
    
    std::cout << "   Debug topo size4: " << topo.nodes.size() << ", edges=" << topo.edges.size() << std::endl;

    // Phase 3: Merge nearby nodes connected by short edges
    //const double merge_threshold = 10.0 * resolution;  // Merge nodes closer than 1 pixel (in meters)
    topo = mergeNearbyNodes(topo, merge_threshold, resolution);

    std::cout << "   Debug topo size5: " << topo.nodes.size() << ", edges=" << topo.edges.size() << std::endl;

    max_search_radius = 100.0 * resolution;
    topo = connectEndpoints(grid, topo, max_search_radius, false);
    
    std::cout << "   Debug topo size6: " << topo.nodes.size() << ", edges=" << topo.edges.size() << std::endl;

    return topo;
}

std::string toJson(const TopologicalMap& map) {
    std::ostringstream os;
    os << "{\n  \"nodes\": [\n";
    for (size_t i = 0; i < map.nodes.size(); ++i) {
        const auto& n = map.nodes[i];
        os << "    {\"id\": " << n.id << ", \"x\": " << n.x << ", \"y\": " << n.y << "}";
        if (i + 1 < map.nodes.size()) os << ",";
        os << "\n";
    }
    os << "  ],\n  \"edges\": [\n";
    for (size_t i = 0; i < map.edges.size(); ++i) {
        const auto& e = map.edges[i];
        os << "    {\"id\": " << e.id << ", \"u\": " << e.u << ", \"v\": " << e.v << ", \"length\": " << e.length ;//<< ", \"polyline\": [";
        /*
        for (size_t j = 0; j < e.polyline.size(); ++j) {
            os << "[" << e.polyline[j].first << ", " << e.polyline[j].second << "]";
            if (j + 1 < e.polyline.size()) os << ", ";
        }
        */
        ///os << "]}";
        os << "}";
        if (i + 1 < map.edges.size()) os << ",";
        os << "\n";
    }
    os << "  ]\n}";
    return os.str();
}

} // namespace gvd_topo


