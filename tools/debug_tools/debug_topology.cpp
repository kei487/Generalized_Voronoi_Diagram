#include <iostream>
#include "gvd_topo/core/OccupancyGrid.hpp"
#include "gvd_topo/core/GvdGenerator.hpp"
#include "gvd_topo/core/TopologyExtractor.hpp"
#include "gvd_topo/io/YamlLoader.hpp"

int main() {
    using namespace gvd_topo;
    
    std::cout << "Debugging topology extraction..." << std::endl;
    
    // Load map
    auto grid = YamlLoader::loadFromYaml("/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/tsudanuma-challenge_nav_resize1.yaml", 50);
    std::cout << "Loaded grid: " << grid.width << "x" << grid.height << std::endl;
    
    // Generate GVD
    GvdGenerator gvd;
    auto result = gvd.run(grid);
    
    size_t gvd_count = 0;
    for (auto v : result.gvd_mask) if (v) ++gvd_count;
    std::cout << "GVD pixels: " << gvd_count << std::endl;
    
    // Count skeleton connectivity
    const int dx8[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
    const int dy8[8] = { -1,-1,-1,  0, 0,  1, 1, 1 };
    auto idx = [](int x, int y, int w) { return y * w + x; };
    auto inBounds = [&](int x, int y) { return x >= 0 && y >= 0 && x < result.width && y < result.height; };
    
    int connectivity_stats[10] = {0}; // 0-9 neighbors
    for (int y = 1; y < result.height-1; ++y) {
        for (int x = 1; x < result.width-1; ++x) {
            if (!result.gvd_mask[idx(x,y,result.width)]) continue;
            int neighbors = 0;
            for (int k = 0; k < 8; ++k) {
                int nx = x + dx8[k];
                int ny = y + dy8[k];
                if (inBounds(nx,ny) && result.gvd_mask[idx(nx,ny,result.width)]) {
                    neighbors++;
                }
            }
            if (neighbors < 10) connectivity_stats[neighbors]++;
        }
    }
    
    std::cout << "Skeleton connectivity stats:" << std::endl;
    for (int i = 0; i < 10; ++i) {
        if (connectivity_stats[i] > 0) {
            std::cout << "  " << i << " neighbors: " << connectivity_stats[i] << " pixels" << std::endl;
        }
    }
    
    // Extract topology with different parameters
    TopologyExtractor::Params params;
    params.prune_min_length = 0.1; // Reduce minimum length
    params.merge_radius = 0.1;     // Reduce merge radius
    
    TopologyExtractor topo(params);
    auto map = topo.run(grid, result.gvd_mask, result.width, result.height, 0.3);
    
    std::cout << "Topological map: " << map.nodes.size() << " nodes, " << map.edges.size() << " edges" << std::endl;
    
    // Analyze node distribution
    if (!map.nodes.empty()) {
        double min_x = map.nodes[0].x, max_x = map.nodes[0].x;
        double min_y = map.nodes[0].y, max_y = map.nodes[0].y;
        for (const auto& node : map.nodes) {
            min_x = std::min(min_x, node.x);
            max_x = std::max(max_x, node.x);
            min_y = std::min(min_y, node.y);
            max_y = std::max(max_y, node.y);
        }
        std::cout << "Node bounds: x[" << min_x << ", " << max_x << "], y[" << min_y << ", " << max_y << "]" << std::endl;
    }
    
    return 0;
}

