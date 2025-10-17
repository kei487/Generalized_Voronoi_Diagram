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

// Compute degree (number of skeleton neighbors) for each pixel
static std::vector<uint8_t> computeDegrees(const std::vector<uint8_t>& is_skel, int width, int height) {
    std::vector<uint8_t> degree(static_cast<size_t>(width * height), 0);
    auto inBounds = [&](int x, int y) { return x >= 0 && y >= 0 && x < width && y < height; };
    
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
    return degree;
}

// Extract raw node candidates (endpoints and junctions)
static std::vector<NodePix> extractRawNodes(const std::vector<uint8_t>& is_skel, const std::vector<uint8_t>& degree, int width, int height) {
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
    
    return raw_nodes;
}

// Merge nearby nodes using spatial grid hashing
static std::vector<int> mergeNearbyNodes(const std::vector<NodePix>& raw_nodes, double merge_radius_px) {
    const double merge_radius_px2 = merge_radius_px * merge_radius_px;
    std::vector<int> parent(raw_nodes.size());
    for (size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
    
    auto findp = [&](int a){ while (parent[a] != a) a = parent[a] = parent[parent[a]]; return a; };
    auto un = [&](int a, int b){ a = findp(a); b = findp(b); if (a!=b) parent[b]=a; };
    
    // Spatial grid optimization
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
        for (int dgy = -1; dgy <= 1; ++dgy) {
            for (int dgx = -1; dgx <= 1; ++dgx) {
                auto it = spatial_grid.find({gx + dgx, gy + dgy});
                if (it == spatial_grid.end()) continue;
                for (int j : it->second) {
                    if (static_cast<size_t>(j) <= i) continue;
                    double dx = static_cast<double>(raw_nodes[i].x - raw_nodes[j].x);
                    double dy = static_cast<double>(raw_nodes[i].y - raw_nodes[j].y);
                    if (dx*dx + dy*dy <= merge_radius_px2) {
                        un(static_cast<int>(i), j);
                    }
                }
            }
        }
    }
    
    return parent;
}

// Create temporary nodes from merged groups
static std::pair<std::vector<TempNode>, std::vector<int>> createTempNodes(
    const std::vector<NodePix>& raw_nodes,
    const std::vector<int>& parent,
    int width, int height,
    double resolution,
    int max_nodes) {
    
    auto findp = [&](int a){ 
        int orig_a = a;
        while (parent[a] != a) a = parent[a];
        return a;
    };
    
    std::vector<std::vector<int>> groups(raw_nodes.size());
    for (size_t i = 0; i < raw_nodes.size(); ++i) {
        groups[findp(static_cast<int>(i))].push_back(static_cast<int>(i));
    }
    
    std::vector<TempNode> temp_nodes;
    temp_nodes.reserve(groups.size());
    std::vector<int> label(static_cast<size_t>(width) * static_cast<size_t>(height), -1);
    
    int temp_id = 0;
    for (const auto& g : groups) if (!g.empty()) {
        if (temp_id >= max_nodes) break;
        double sx=0, sy=0;
        for (int id : g) { sx += raw_nodes[id].x; sy += raw_nodes[id].y; }
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
    
    return {temp_nodes, label};
}

// Filter nodes by connected component analysis on skeleton
static std::pair<std::vector<TopoNode>, std::vector<int>> filterByConnectedComponents(
    const std::vector<TempNode>& temp_nodes,
    const std::vector<uint8_t>& is_skel,
    std::vector<int> label,
    int width, int height,
    int min_component_size) {
    
    auto inBounds = [&](int x, int y) { return x >= 0 && y >= 0 && x < width && y < height; };
    
    std::vector<bool> node_visited(temp_nodes.size(), false);
    std::vector<int> component_id(temp_nodes.size(), -1);
    std::vector<int> component_size;
    int num_components = 0;
    
    for (size_t start_idx = 0; start_idx < temp_nodes.size(); ++start_idx) {
        if (node_visited[start_idx]) continue;
        
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
            
            std::unordered_set<int> local_visited_set;
            std::queue<std::pair<int,int>> skel_q;
            skel_q.push({nx, ny});
            local_visited_set.insert(idx(nx, ny, width));
            const int max_search_dist = 50;
            int search_count = 0;
            
            while (!skel_q.empty() && search_count < max_search_dist) {
                auto [cx, cy] = skel_q.front();
                skel_q.pop();
                ++search_count;
                
                for (int k = 0; k < 8; ++k) {
                    int tx = cx + dx8[k];
                    int ty = cy + dy8[k];
                    if (!inBounds(tx, ty) || !is_skel[idx(tx, ty, width)]) continue;
                    int tidx = idx(tx, ty, width);
                    if (local_visited_set.count(tidx)) continue;
                    local_visited_set.insert(tidx);
                    
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
    
    // Find largest component and filter
    int largest_comp = 0;
    int largest_size = component_size.empty() ? 0 : component_size[0];
    for (int i = 1; i < num_components; ++i) {
        if (component_size[i] > largest_size) {
            largest_size = component_size[i];
            largest_comp = i;
        }
    }
    
    std::vector<bool> keep_component(num_components, false);
    for (int i = 0; i < num_components; ++i) {
        if (component_size[i] > min_component_size || i == largest_comp) {
            keep_component[i] = true;
        }
    }
    
    // Rebuild nodes and label map
    std::vector<TopoNode> nodes;
    label.assign(static_cast<size_t>(width) * static_cast<size_t>(height), -1);
    int node_id = 0;
    
    for (size_t i = 0; i < temp_nodes.size(); ++i) {
        int comp_id = component_id[i];
        if (comp_id >= 0 && comp_id < num_components && keep_component[comp_id]) {
            TopoNode n;
            n.id = node_id;
            n.x = temp_nodes[i].world_x;
            n.y = temp_nodes[i].world_y;
            nodes.push_back(n);
            label[idx(temp_nodes[i].x, temp_nodes[i].y, width)] = node_id;
            ++node_id;
        }
    }
    
    return {nodes, label};
}

// Trace edges from nodes along skeleton
struct EdgeTracingResult {
    std::vector<TopoEdge> edges;
    std::vector<bool> node_used;
    int new_nodes_created;
};

static EdgeTracingResult traceEdges(
    std::vector<TopoNode>& nodes,
    const std::vector<uint8_t>& is_skel,
    std::vector<int>& label,
    int width, int height,
    double resolution,
    int max_edges) {
    
    auto inBounds = [&](int x, int y) { return x >= 0 && y >= 0 && x < width && y < height; };
    
    std::vector<uint8_t> visited(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
    std::vector<bool> node_used(nodes.size(), false);
    std::vector<TopoEdge> edges;
    
    auto isNode = [&](int x, int y) {
        if (!inBounds(x, y)) return false;
        int id = label[idx(x,y,width)];
        return id >= 0;
    };
    
    int edge_count = 0;
    int node_id = static_cast<int>(nodes.size());
    const int max_nodes = 100000;
    
    for (const auto& n : nodes) {
        if (edge_count >= max_edges) break;
        int sx = static_cast<int>(std::round(n.x / resolution));
        int sy = static_cast<int>(std::round(n.y / resolution));
        
        // Find nearest skeleton pixel if not on skeleton
        if (!inBounds(sx, sy) || !is_skel[idx(sx, sy, width)]) {
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
        
        // Trace in 8 directions
        for (int k = 0; k < 8; ++k) {
            int nx = sx + dx8[k];
            int ny = sy + dy8[k];
            if (!inBounds(nx,ny) || !is_skel[idx(nx,ny,width)]) continue;
            if (visited[idx(nx,ny,width)]) continue;
            
            // Trace path
            std::vector<std::pair<double,double>> poly;
            int px = sx; int py = sy; int cx = nx; int cy = ny;
            double length = 0.0;
            int steps = 0;
            const int max_steps = std::min(width * height, 100000);
            
            while (steps < max_steps) {
                visited[idx(cx,cy,width)] = 1;
                poly.emplace_back(cx * resolution, cy * resolution);
                
                if (isNode(cx,cy) && !(cx == sx && cy == sy)) {
                    int to_id = label[idx(cx,cy,width)];
                    TopoEdge e;
                    e.id = static_cast<int>(edges.size());
                    e.u = n.id;
                    e.v = to_id;
                    e.length = length;
                    e.polyline = poly;
                    edges.push_back(std::move(e));
                    // Mark both nodes as used
                    if (n.id < static_cast<int>(node_used.size())) node_used[n.id] = true;
                    if (to_id < static_cast<int>(node_used.size())) node_used[to_id] = true;
                    ++edge_count;
                    break;
                }
                
                if (cx == sx && cy == sy) break;
                
                // Choose next neighbor
                int nextx = -1, nexty = -1, choices = 0;
                for (int kk = 0; kk < 8; ++kk) {
                    int tx = cx + dx8[kk];
                    int ty = cy + dy8[kk];
                    if (!inBounds(tx,ty) || !is_skel[idx(tx,ty,width)]) continue;
                    if (tx == px && ty == py) continue;
                    if (visited[idx(tx,ty,width)]) continue;
                    ++choices; nextx = tx; nexty = ty;
                }
                
                if (choices == 0) {
                    // Dead end -> create node if not exist
                    if (!isNode(cx,cy) && node_id < max_nodes) {
                        TopoNode m;
                        m.id = node_id;
                        m.x = cx * resolution;
                        m.y = cy * resolution;
                        nodes.push_back(m);
                        node_used.push_back(true);
                        label[idx(cx,cy,width)] = node_id;
                        ++node_id;
                        
                        int to_id = label[idx(cx,cy,width)];
                        TopoEdge e;
                        e.id = static_cast<int>(edges.size());
                        e.u = n.id;
                        e.v = to_id;
                        e.length = length;
                        e.polyline = poly;
                        edges.push_back(std::move(e));
                        if (n.id < static_cast<int>(node_used.size())) node_used[n.id] = true;
                        ++edge_count;
                    }
                    break;
                }
                
                length += std::hypot(static_cast<double>(nextx - cx) * resolution,
                                    static_cast<double>(nexty - cy) * resolution);
                px = cx; py = cy; cx = nextx; cy = nexty;
                ++steps;
            }
        }
    }
    
    return {edges, node_used, node_id - static_cast<int>(nodes.size())};
}

// Remove unused nodes after edge tracing
static std::pair<std::vector<TopoNode>, std::vector<TopoEdge>> removeUnusedNodes(
    const std::vector<TopoNode>& nodes,
    const std::vector<TopoEdge>& edges,
    const std::vector<bool>& node_used) {
    
    std::vector<TopoNode> used_nodes;
    used_nodes.reserve(nodes.size());
    std::map<int, int> old_to_new_node_id;
    int new_node_id = 0;
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (i < node_used.size() && node_used[i]) {
            TopoNode new_node = nodes[i];
            new_node.id = new_node_id;
            old_to_new_node_id[nodes[i].id] = new_node_id;
            used_nodes.push_back(new_node);
            ++new_node_id;
        }
    }
    
    // Update edge node IDs
    std::vector<TopoEdge> updated_edges;
    updated_edges.reserve(edges.size());
    
    for (auto e : edges) {
        if (old_to_new_node_id.find(e.u) != old_to_new_node_id.end() &&
            old_to_new_node_id.find(e.v) != old_to_new_node_id.end()) {
            e.u = old_to_new_node_id[e.u];
            e.v = old_to_new_node_id[e.v];
            e.id = static_cast<int>(updated_edges.size());
            updated_edges.push_back(e);
        }
    }
    
    return {used_nodes, updated_edges};
}

// Prune short edges (dead-end spurs)
static std::vector<TopoEdge> pruneShortEdges(
    const std::vector<TopoEdge>& edges,
    double min_length) {
    
    // Build degree map
    std::map<int, int> node_degree;
    for (const auto& e : edges) {
        node_degree[e.u]++;
        node_degree[e.v]++;
    }
    
    std::vector<TopoEdge> kept;
    kept.reserve(edges.size());
    
    for (const auto& e : edges) {
        bool should_keep = true;
        if (e.length < min_length) {
            bool u_is_endpoint = (node_degree[e.u] <= 1);
            bool v_is_endpoint = (node_degree[e.v] <= 1);
            should_keep = !(u_is_endpoint && v_is_endpoint);
        }
        if (should_keep) kept.push_back(e);
    }
    
    return kept;
}

// ============================================================================
// Main Topology Extraction Function
// ============================================================================

TopologicalMap TopologyExtractor::run(const std::vector<uint8_t>& gvd_mask, int width, int height, double resolution) const {
    TopologicalMap topo;
    if (gvd_mask.empty() || width <= 0 || height <= 0) return topo;

    // Phase 0: Apply Zhang-Suen thinning as preprocessing
    std::vector<uint8_t> is_skel = zhangSuenThinning(gvd_mask, width, height);

    // Phase 1: Compute degrees for each skeleton pixel
    std::vector<uint8_t> degree = computeDegrees(is_skel, width, height);

    // Phase 2: Extract raw node candidates (endpoints and junctions)
    std::vector<NodePix> raw_nodes = extractRawNodes(is_skel, degree, width, height);

    // Phase 3: Merge nearby nodes using spatial hashing
    const double merge_radius_px = params_.merge_radius / resolution;
    std::vector<int> parent = mergeNearbyNodes(raw_nodes, merge_radius_px);

    // Phase 4: Create temporary nodes from merged groups
    const int max_nodes = 100000;
    auto [temp_nodes, label] = createTempNodes(raw_nodes, parent, width, height, resolution, max_nodes);

    // Phase 5: Filter by connected components
    /*
    const size_t max_nodes_for_component_analysis = 10000;
    const int min_component_size = 1;
    
    auto [filtered_nodes, updated_label] = filterByConnectedComponents(
        temp_nodes, is_skel, label, width, height, min_component_size);
    topo.nodes = filtered_nodes;
    label = updated_label;
    */
    
    // Phase 6: Trace edges along skeleton
    const int max_edges = 200000;
    auto [traced_edges, node_used, new_nodes_count] = traceEdges(
        topo.nodes, is_skel, label, width, height, resolution, max_edges);
    topo.edges = traced_edges;

    // Phase 7: Remove unused nodes after edge tracing
    auto [used_nodes, updated_edges] = removeUnusedNodes(topo.nodes, topo.edges, node_used);
    topo.nodes = used_nodes;
    topo.edges = updated_edges;

    // Phase 8: Prune short dead-end edges
    topo.edges = pruneShortEdges(topo.edges, params_.prune_min_length);

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


