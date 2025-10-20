#include <iostream>
#include "gvd_topo/core/OccupancyGrid.hpp"
#include "gvd_topo/core/GvdGenerator.hpp"
#include "gvd_topo/core/TopologyExtractor.hpp"
#include "gvd_topo/io/YamlLoader.hpp"
#include "gvd_topo/core/Visualizer.hpp"

int main() {
    using namespace gvd_topo;
    
    std::cout << "Debugging visualization..." << std::endl;
    
    // Load map
    auto grid = YamlLoader::loadFromYaml("/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/tsudanuma-challenge_nav_resize1.yaml", 50);
    std::cout << "Loaded grid: " << grid.width << "x" << grid.height << ", resolution: " << grid.resolution << std::endl;
    std::cout << "Origin: (" << grid.origin.x << ", " << grid.origin.y << ", " << grid.origin.theta << ")" << std::endl;
    
    // Generate GVD
    GvdGenerator gvd;
    auto result = gvd.run(grid);
    
    size_t gvd_count = 0;
    for (auto v : result.gvd_mask) if (v) ++gvd_count;
    std::cout << "GVD pixels: " << gvd_count << std::endl;
    
    // Extract topology
    TopologyExtractor topo;
    auto map = topo.run(grid, result.gvd_mask, result.width, result.height, grid.resolution);
    
    std::cout << "Using resolution: " << grid.resolution << " m/pixel" << std::endl;
    
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
    
    // Analyze edge distribution
    if (!map.edges.empty()) {
        double min_len = map.edges[0].length;
        double max_len = map.edges[0].length;
        int edges_with_polyline = 0;
        for (const auto& edge : map.edges) {
            min_len = std::min(min_len, edge.length);
            max_len = std::max(max_len, edge.length);
            //if (!edge.polyline.empty()) edges_with_polyline++;
        }
        std::cout << "Edge length range: [" << min_len << ", " << max_len << "]" << std::endl;
        std::cout << "Edges with polyline: " << edges_with_polyline << "/" << map.edges.size() << std::endl;
    }
    
    // Test visualization with different options
    VisualizationOptions options;
    options.image_width = 1200;
    options.image_height = 800;
    options.node_radius = 5;
    options.edge_thickness = 2;
    options.show_node_ids = true;
    options.show_edge_lengths = true;
    
    try {
        Visualizer::saveTopologicalMapAsImage(map, "/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/debug_topo.png", options);
        std::cout << "Debug visualization saved to debug_topo.png" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Visualization error: " << e.what() << std::endl;
    }
    
    return 0;
}
