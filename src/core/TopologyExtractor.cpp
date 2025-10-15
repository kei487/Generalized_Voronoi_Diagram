#include "gvd_topo/core/TopologyExtractor.hpp"
#include <cstdint>
#include <queue>
#include <cmath>
#include <limits>
#include <sstream>
#include <map>
#include <algorithm>
#include <vector>
#include <unordered_set>

namespace gvd_topo {

TopologyExtractor::TopologyExtractor() = default;
TopologyExtractor::TopologyExtractor(const Params& p) : params_(p) {}

static inline int idx(int x, int y, int w) { return y * w + x; }

// Optimized Zhang-Suen thinning algorithm with parallel processing
static std::vector<uint8_t> zhangSuenThinning(const std::vector<uint8_t>& input, int width, int height) {
    std::vector<uint8_t> img = input;
    std::vector<uint8_t> marker(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
    
    bool changed = true;
    int iteration = 0;
    // Adjust max iterations based on map size
    const int max_iterations = (width * height > 10000000) ? 20 : 100;
    
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

TopologicalMap TopologyExtractor::run(const std::vector<uint8_t>& gvd_mask, int width, int height, double resolution) const {
    TopologicalMap topo;
    if (gvd_mask.empty() || width <= 0 || height <= 0) return topo;

    auto inBounds = [&](int x, int y) { return x >= 0 && y >= 0 && x < width && y < height; };

    // Apply Zhang-Suen thinning as preprocessing
    // Skip thinning for very large maps to avoid excessive processing time
    const int max_pixels_for_thinning = 5000000; // 5 million pixels
    std::vector<uint8_t> is_skel;
    is_skel = zhangSuenThinning(gvd_mask, width, height);



    // Compute degrees for each skeleton pixel (8-neighborhood)
    std::vector<uint8_t> degree(static_cast<size_t>(width * height), 0);
    const int dx8[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
    const int dy8[8] = { -1,-1,-1,  0, 0,  1, 1, 1 };
    #ifdef GVD_TOPO_WITH_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (!is_skel[idx(x,y,width)]) continue;
            int deg = 0;
            for (int k = 0; k < 8; ++k) {
                int nx = x + dx8[k];
                int ny = y + dy8[k];
                if (inBounds(nx,ny) && is_skel[idx(nx,ny,width)]) ++deg;
            }
            degree[idx(x,y,width)] = static_cast<uint8_t>(deg);
        }
    }

    // Identify raw nodes (endpoints degree==1, junctions degree>=3)
    struct NodePix { int x; int y; };
    std::vector<NodePix> raw_nodes;
    raw_nodes.reserve(1024);
    #ifdef GVD_TOPO_WITH_OPENMP
    #pragma omp parallel
    {
        std::vector<NodePix> local_nodes;
        local_nodes.reserve(1024);
        #pragma omp for nowait schedule(static)
        for (int y = 1; y < height-1; ++y) {
            for (int x = 1; x < width-1; ++x) {
                if (!is_skel[idx(x,y,width)]) continue;
                int d = degree[idx(x,y,width)];
                if (d == 1 || d >= 3) local_nodes.push_back({x,y});
            }
        }
        #pragma omp critical
        raw_nodes.insert(raw_nodes.end(), local_nodes.begin(), local_nodes.end());
    }
    #else
    for (int y = 1; y < height-1; ++y) {
        for (int x = 1; x < width-1; ++x) {
            if (!is_skel[idx(x,y,width)]) continue;
            int d = degree[idx(x,y,width)];
            if (d == 1 || d >= 3) raw_nodes.push_back({x,y});
        }
    }
    #endif

    // Merge nearby nodes within merge_radius (pixels)
    // Use grid-based spatial hashing to avoid O(N²) comparison
    const double merge_radius_px = params_.merge_radius / resolution;
    const double merge_radius_px2 = merge_radius_px * merge_radius_px;
    std::vector<int> parent(raw_nodes.size());
    for (size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
    auto findp = [&](int a){ while (parent[a] != a) a = parent[a] = parent[parent[a]]; return a; };
    auto un = [&](int a, int b){ a = findp(a); b = findp(b); if (a!=b) parent[b]=a; };
    
    // Spatial grid optimization: only compare nodes in nearby grid cells
    const int grid_size = static_cast<int>(std::max(1.0, merge_radius_px)) + 1;
    std::map<std::pair<int,int>, std::vector<int>> spatial_grid;
    for (size_t i = 0; i < raw_nodes.size(); ++i) {
        int gx = raw_nodes[i].x / grid_size;
        int gy = raw_nodes[i].y / grid_size;
        spatial_grid[{gx, gy}].push_back(static_cast<int>(i));
    }
    
    // Compare only nodes in same or adjacent grid cells
    for (size_t i = 0; i < raw_nodes.size(); ++i) {
        int gx = raw_nodes[i].x / grid_size;
        int gy = raw_nodes[i].y / grid_size;
        // Check 3x3 neighborhood of grid cells
        for (int dgy = -1; dgy <= 1; ++dgy) {
            for (int dgx = -1; dgx <= 1; ++dgx) {
                auto it = spatial_grid.find({gx + dgx, gy + dgy});
                if (it == spatial_grid.end()) continue;
                for (int j : it->second) {
                    if (static_cast<size_t>(j) <= i) continue; // Avoid duplicate comparisons
                    double dx = static_cast<double>(raw_nodes[i].x - raw_nodes[j].x);
                    double dy = static_cast<double>(raw_nodes[i].y - raw_nodes[j].y);
                    if (dx*dx + dy*dy <= merge_radius_px2) {
                        un(static_cast<int>(i), j);
                    }
                }
            }
        }
    }
    // Compute representatives and average positions
    std::vector<std::vector<int>> groups(raw_nodes.size());
    for (size_t i = 0; i < raw_nodes.size(); ++i) groups[findp(static_cast<int>(i))].push_back(static_cast<int>(i));
    
    // First pass: create temporary nodes and check connectivity
    struct TempNode {
        int id;
        int x, y;  // pixel coordinates
        double world_x, world_y;
    };
    std::vector<TempNode> temp_nodes;
    temp_nodes.reserve(groups.size());
    std::vector<int> label(static_cast<size_t>(width) * static_cast<size_t>(height), -1);
    
    const int max_nodes = 100000;
    int temp_id = 0;
    for (const auto& g : groups) if (!g.empty()) {
        if (temp_id >= max_nodes) break;
        double sx=0, sy=0; 
        for (int id : g){ sx += raw_nodes[id].x; sy += raw_nodes[id].y; }
        int cx = static_cast<int>(std::round(sx / static_cast<double>(g.size())));
        int cy = static_cast<int>(std::round(sy / static_cast<double>(g.size())));
        if (cx < 0 || cx >= width || cy < 0 || cy >= height) continue;
        
        TempNode tn;
        tn.id = temp_id;
        tn.x = cx;
        tn.y = cy;
        tn.world_x = cx * resolution;
        tn.world_y = cy * resolution;
        temp_nodes.push_back(tn);
        label[idx(cx,cy,width)] = temp_id;
        ++temp_id;
    }
    
    // Connected component analysis: use simplified method for large maps
    // For large maps (>10K nodes), skip this step and rely on edge tracing to filter nodes
    const size_t max_nodes_for_component_analysis = 10000;
    
    if (temp_nodes.size() <= max_nodes_for_component_analysis) {
        // Use BFS to find connected components (for smaller maps only)
        std::vector<bool> node_visited(temp_nodes.size(), false);
        std::vector<int> component_id(temp_nodes.size(), -1);
        std::vector<int> component_size;
        int num_components = 0;
        
        for (size_t start_idx = 0; start_idx < temp_nodes.size(); ++start_idx) {
            if (node_visited[start_idx]) continue;
            
            // BFS from this node along skeleton
            std::queue<int> q;
            q.push(static_cast<int>(start_idx));
            node_visited[start_idx] = true;
            component_id[start_idx] = num_components;
            int comp_size = 0;
            
            while (!q.empty()) {
                int curr_node_idx = q.front();
                q.pop();
                ++comp_size;
                
                int nx = temp_nodes[curr_node_idx].x;
                int ny = temp_nodes[curr_node_idx].y;
                
                // Search along skeleton to find neighboring nodes (limited search)
                std::unordered_set<int> local_visited_set;
                std::queue<std::pair<int,int>> skel_q;
                skel_q.push({nx, ny});
                local_visited_set.insert(idx(nx, ny, width));
                const int max_search_dist = 50; // Reduced for performance
                int search_count = 0;
                
                while (!skel_q.empty() && search_count < max_search_dist) {
                    auto [cx, cy] = skel_q.front();
                    skel_q.pop();
                    ++search_count;
                    
                    // Check 8 neighbors on skeleton
                    for (int k = 0; k < 8; ++k) {
                        int tx = cx + dx8[k];
                        int ty = cy + dy8[k];
                        if (!inBounds(tx, ty) || !is_skel[idx(tx, ty, width)]) continue;
                        int tidx = idx(tx, ty, width);
                        if (local_visited_set.count(tidx)) continue;
                        local_visited_set.insert(tidx);
                        
                        // Check if this position has a node
                        int neighbor_label = label[tidx];
                        if (neighbor_label >= 0 && neighbor_label != temp_nodes[curr_node_idx].id) {
                            size_t neighbor_idx = static_cast<size_t>(neighbor_label);
                            if (neighbor_idx < temp_nodes.size() && !node_visited[neighbor_idx]) {
                                node_visited[neighbor_idx] = true;
                                component_id[neighbor_idx] = num_components;
                                q.push(static_cast<int>(neighbor_idx));
                            }
                        }
                        
                        skel_q.push({tx, ty});
                    }
                }
            }
            
            component_size.push_back(comp_size);
            ++num_components;
        }
        
        // Filter out small connected components (size <= 20) and keep larger ones
        // Always keep the largest component even if it's small
        const int min_component_size = 1;
        if (!component_size.empty()) {
            // Find the largest component
            int largest_comp = 0;
            int largest_size = component_size[0];
            for (int i = 1; i < num_components; ++i) {
                if (component_size[i] > largest_size) {
                    largest_size = component_size[i];
                    largest_comp = i;
                }
            }
            
            // Identify which components to keep (size > threshold OR is the largest)
            std::vector<bool> keep_component(num_components, false);
            for (int i = 0; i < num_components; ++i) {
                if (component_size[i] > min_component_size || i == largest_comp) {
                    keep_component[i] = true;
                }
            }
            
            // Rebuild nodes and label map with only nodes from large enough components
            topo.nodes.clear();
            label.assign(static_cast<size_t>(width) * static_cast<size_t>(height), -1);
            int node_id = 0;
            
            for (size_t i = 0; i < temp_nodes.size(); ++i) {
                int comp_id = component_id[i];
                if (comp_id >= 0 && comp_id < num_components && keep_component[comp_id]) {
                    TopoNode n;
                    n.id = node_id;
                    n.x = temp_nodes[i].world_x;
                    n.y = temp_nodes[i].world_y;
                    topo.nodes.push_back(n);
                    label[idx(temp_nodes[i].x, temp_nodes[i].y, width)] = node_id;
                    ++node_id;
                }
            }
        }
    } else {
        // For large maps, directly use temp_nodes without component analysis
        // Edge tracing will filter out isolated nodes later
        for (const auto& tn : temp_nodes) {
            TopoNode n;
            n.id = topo.nodes.size();
            n.x = tn.world_x;
            n.y = tn.world_y;
            topo.nodes.push_back(n);
        }
    }

    // Edge tracing: from each node, follow skeleton until another node or endpoint
    std::vector<uint8_t> visited(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
    std::vector<bool> node_used(topo.nodes.size(), false); // Track which nodes are actually used
    auto isNode = [&](int x, int y){ 
        if (!inBounds(x, y)) return false;
        int id = label[idx(x,y,width)]; 
        return id >= 0; 
    };

    const int max_edges = 200000; // Limit number of edges to prevent memory issues
    int edge_count = 0;
    int node_id = static_cast<int>(topo.nodes.size()); // Continue node ID from existing nodes
    for (const auto& n : topo.nodes) {
        if (edge_count >= max_edges) break; // Safety limit
        int sx = static_cast<int>(std::round(n.x / resolution));
        int sy = static_cast<int>(std::round(n.y / resolution));
        
        // Check if start position is on skeleton
        if (!inBounds(sx, sy) || !is_skel[idx(sx, sy, width)]) {
            // Find nearest skeleton pixel
            int best_x = sx, best_y = sy, best_dist = 1000;
            for (int dy = -3; dy <= 3; ++dy) {
                for (int dx = -3; dx <= 3; ++dx) {
                    int nx = sx + dx, ny = sy + dy;
                    if (inBounds(nx, ny) && is_skel[idx(nx, ny, width)]) {
                        int dist = dx*dx + dy*dy;
                        if (dist < best_dist) {
                            best_dist = dist;
                            best_x = nx; best_y = ny;
                        }
                    }
                }
            }
            sx = best_x; sy = best_y;
        }
        
        for (int k = 0; k < 8; ++k) {
            int nx = sx + dx8[k];
            int ny = sy + dy8[k];
            if (!inBounds(nx,ny) || !is_skel[idx(nx,ny,width)]) continue;
            if (visited[idx(nx,ny,width)]) continue;
            
            // trace path
            std::vector<std::pair<double,double>> poly;
            int px = sx; int py = sy; int cx = nx; int cy = ny;
            double length = 0.0;
            int steps = 0;
            const int max_steps = std::min(width * height, 100000); // Limit to prevent infinite loops
            while (steps < max_steps) {
                visited[idx(cx,cy,width)] = 1;
                poly.emplace_back(cx * resolution, cy * resolution);
                
                if (isNode(cx,cy) && !(cx == sx && cy == sy)) {
                    int to_id = label[idx(cx,cy,width)];
                    TopoEdge e; e.id = static_cast<int>(topo.edges.size()); e.u = n.id; e.v = to_id; e.length = length; e.polyline = poly;
                    topo.edges.push_back(std::move(e));
                    // Mark both nodes as used
                    if (n.id < static_cast<int>(node_used.size())) node_used[n.id] = true;
                    if (to_id < static_cast<int>(node_used.size())) node_used[to_id] = true;
                    ++edge_count;
                    break;
                }
                
                // prevent cycling back to start node
                if (cx == sx && cy == sy) {
                    break;
                }
                
                // choose next neighbor (avoid returning to previous pixel)
                int nextx = -1, nexty = -1, choices = 0;
                for (int kk = 0; kk < 8; ++kk) {
                    int tx = cx + dx8[kk];
                    int ty = cy + dy8[kk];
                    if (!inBounds(tx,ty) || !is_skel[idx(tx,ty,width)]) continue;
                    if (tx == px && ty == py) continue;
                    if (visited[idx(tx,ty,width)]) continue; // avoid cycles
                    ++choices; nextx = tx; nexty = ty;
                }
                
                if (choices == 0) {
                    // dead end -> create node if not exist
                    if (!isNode(cx,cy) && node_id < max_nodes) {
                        TopoNode m; m.id = node_id; m.x = cx * resolution; m.y = cy * resolution; 
                        topo.nodes.push_back(m); 
                        node_used.push_back(true); // New node is immediately used
                        label[idx(cx,cy,width)] = node_id; 
                        ++node_id;
                        int to_id = label[idx(cx,cy,width)];
                        TopoEdge e; e.id = static_cast<int>(topo.edges.size()); e.u = n.id; e.v = to_id; e.length = length; e.polyline = poly;
                        topo.edges.push_back(std::move(e));
                        // Mark start node as used
                        if (n.id < static_cast<int>(node_used.size())) node_used[n.id] = true;
                        ++edge_count;
                    }
                    break;
                }
                
                length += std::hypot(static_cast<double>(nextx - cx) * resolution, static_cast<double>(nexty - cy) * resolution);
                px = cx; py = cy; cx = nextx; cy = nexty;
                ++steps;
            }
        }
    }

    // Remove unused nodes after edge tracing (second method: actual usage tracking)
    std::vector<TopoNode> used_nodes;
    used_nodes.reserve(topo.nodes.size());
    std::map<int, int> old_to_new_node_id;
    int new_node_id = 0;
    
    for (size_t i = 0; i < topo.nodes.size(); ++i) {
        if (i < node_used.size() && node_used[i]) {
            TopoNode new_node = topo.nodes[i];
            new_node.id = new_node_id;
            old_to_new_node_id[topo.nodes[i].id] = new_node_id;
            used_nodes.push_back(new_node);
            ++new_node_id;
        }
    }
    
    topo.nodes.swap(used_nodes);
    
    // Update edge node IDs
    for (auto& e : topo.edges) {
        if (old_to_new_node_id.find(e.u) != old_to_new_node_id.end()) {
            e.u = old_to_new_node_id[e.u];
        }
        if (old_to_new_node_id.find(e.v) != old_to_new_node_id.end()) {
            e.v = old_to_new_node_id[e.v];
        }
    }
    
    // Additional filtering: remove small connected components after edge tracing
    // Build adjacency list from edges
    std::map<int, std::vector<int>> adjacency;
    for (const auto& e : topo.edges) {
        adjacency[e.u].push_back(e.v);
        adjacency[e.v].push_back(e.u);
    }
    
    // Find connected components using BFS on the graph
    std::vector<bool> visited_nodes(topo.nodes.size(), false);
    std::vector<int> node_component(topo.nodes.size(), -1);
    std::vector<int> comp_sizes;
    int comp_count = 0;
    
    for (size_t i = 0; i < topo.nodes.size(); ++i) {
        if (visited_nodes[i]) continue;
        
        // BFS to find all nodes in this component
        std::queue<int> q;
        q.push(static_cast<int>(i));
        visited_nodes[i] = true;
        node_component[i] = comp_count;
        int size = 0;
        
        while (!q.empty()) {
            int curr = q.front();
            q.pop();
            ++size;
            
            if (adjacency.find(curr) != adjacency.end()) {
                for (int neighbor : adjacency[curr]) {
                    if (neighbor >= 0 && neighbor < static_cast<int>(topo.nodes.size()) && !visited_nodes[neighbor]) {
                        visited_nodes[neighbor] = true;
                        node_component[neighbor] = comp_count;
                        q.push(neighbor);
                    }
                }
            }
        }
        
        comp_sizes.push_back(size);
        ++comp_count;
    }
    
    // Filter out small components (size <= 20), but always keep the largest
    const int min_final_component_size = 20;
    
    // Find the largest component
    int largest_final_comp = 0;
    int largest_final_size = comp_count > 0 ? comp_sizes[0] : 0;
    for (int i = 1; i < comp_count; ++i) {
        if (comp_sizes[i] > largest_final_size) {
            largest_final_size = comp_sizes[i];
            largest_final_comp = i;
        }
    }
    
    std::vector<bool> keep_final_component(comp_count, false);
    for (int i = 0; i < comp_count; ++i) {
        if (comp_sizes[i] > min_final_component_size || i == largest_final_comp) {
            keep_final_component[i] = true;
        }
    }
    
    // Rebuild nodes and edges, keeping only nodes from large components
    std::vector<TopoNode> final_nodes;
    final_nodes.reserve(topo.nodes.size());
    std::map<int, int> final_old_to_new;
    int final_node_id = 0;
    
    for (size_t i = 0; i < topo.nodes.size(); ++i) {
        int comp_id = node_component[i];
        if (comp_id >= 0 && comp_id < comp_count && keep_final_component[comp_id]) {
            TopoNode n = topo.nodes[i];
            n.id = final_node_id;
            final_old_to_new[static_cast<int>(i)] = final_node_id;
            final_nodes.push_back(n);
            ++final_node_id;
        }
    }
    
    topo.nodes.swap(final_nodes);
    
    // Update edges with new node IDs and remove edges connecting to removed nodes
    std::vector<TopoEdge> final_edges;
    final_edges.reserve(topo.edges.size());
    
    for (auto& e : topo.edges) {
        if (final_old_to_new.find(e.u) != final_old_to_new.end() &&
            final_old_to_new.find(e.v) != final_old_to_new.end()) {
            e.u = final_old_to_new[e.u];
            e.v = final_old_to_new[e.v];
            e.id = static_cast<int>(final_edges.size());
            final_edges.push_back(e);
        }
    }
    
    topo.edges.swap(final_edges);

    // Pruning: remove edges shorter than threshold
    const double min_len = params_.prune_min_length;
    
    // Build degree map efficiently (O(E) instead of O(E²))
    std::map<int, int> node_degree;
    for (const auto& e : topo.edges) {
        node_degree[e.u]++;
        node_degree[e.v]++;
    }
    
    std::vector<TopoEdge> kept;
    kept.reserve(topo.edges.size());
    for (auto& e : topo.edges) {
        // Only prune if edge is very short AND connects to degree-1 nodes (endpoints)
        bool should_keep = true;
        if (e.length < min_len) {
            // Check if both endpoints are degree-1 nodes (dead ends)
            bool u_is_endpoint = (node_degree[e.u] <= 1);
            bool v_is_endpoint = (node_degree[e.v] <= 1);
            // Only prune if both endpoints are dead ends
            should_keep = !(u_is_endpoint && v_is_endpoint);
        }
        if (should_keep) kept.push_back(std::move(e));
    }
    topo.edges.swap(kept);

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
        os << "    {\"id\": " << e.id << ", \"u\": " << e.u << ", \"v\": " << e.v << ", \"length\": " << e.length << ", \"polyline\": [";
        for (size_t j = 0; j < e.polyline.size(); ++j) {
            os << "[" << e.polyline[j].first << ", " << e.polyline[j].second << "]";
            if (j + 1 < e.polyline.size()) os << ", ";
        }
        os << "]}";
        if (i + 1 < map.edges.size()) os << ",";
        os << "\n";
    }
    os << "  ]\n}";
    return os.str();
}

} // namespace gvd_topo


