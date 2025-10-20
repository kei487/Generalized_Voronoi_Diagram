#include <iostream>
#include "gvd_topo/core/OccupancyGrid.hpp"
#include "gvd_topo/core/GvdGenerator.hpp"
#include "gvd_topo/core/TopologyExtractor.hpp"
#include "gvd_topo/io/YamlLoader.hpp"
#include "gvd_topo/core/Visualizer.hpp"

int main() {
    using namespace gvd_topo;
    
    std::cout << "Debugging new map..." << std::endl;
    
    // Load map
    auto grid = YamlLoader::loadFromYaml("/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/map_tsudanuma.yaml", 50);
    std::cout << "Loaded grid: " << grid.width << "x" << grid.height << ", resolution: " << grid.resolution << std::endl;
    std::cout << "Origin: (" << grid.origin.x << ", " << grid.origin.y << ", " << grid.origin.theta << ")" << std::endl;
    
    // Count occupied pixels
    int occupied = 0;
    for (auto v : grid.data) if (v > 50) ++occupied;
    std::cout << "Occupied pixels: " << occupied << "/" << (grid.width * grid.height) << std::endl;
    
    // Generate GVD
    std::cout << "Generating GVD..." << std::endl;
    GvdGenerator gvd;
    auto result = gvd.run(grid);
    
    size_t gvd_count = 0;
    for (auto v : result.gvd_mask) if (v) ++gvd_count;
    std::cout << "GVD pixels: " << gvd_count << std::endl;
    
    // Extract topology
    std::cout << "Extracting topology..." << std::endl;
    TopologyExtractor topo;
    auto map = topo.run(grid, result.gvd_mask, result.width, result.height, grid.resolution);
    
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
        std::cout << "Node range: " << (max_x - min_x) << " x " << (max_y - min_y) << std::endl;
    }
    
    // Test visualization with smaller image
    VisualizationOptions options;
    options.image_width = 800;
    options.image_height = 600;
    options.node_radius = 3;
    options.edge_thickness = 1;
    options.show_node_ids = false;
    options.show_edge_lengths = false;
    
    try {
        Visualizer::saveTopologicalMapAsImage(map, "/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/debug_new_topo.png", options);
        std::cout << "Debug visualization saved to debug_new_topo.png" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Visualization error: " << e.what() << std::endl;
    }
    
    return 0;
}

