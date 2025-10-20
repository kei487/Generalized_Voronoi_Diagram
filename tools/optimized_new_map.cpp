#include <iostream>
#include <fstream>
#include "gvd_topo/core/OccupancyGrid.hpp"
#include "gvd_topo/core/GvdGenerator.hpp"
#include "gvd_topo/core/TopologyExtractor.hpp"
#include "gvd_topo/io/YamlLoader.hpp"
#include "gvd_topo/core/Visualizer.hpp"

int main() {
    using namespace gvd_topo;
    
    std::cout << "Processing large map with optimized parameters..." << std::endl;
    
    // Load map
    auto grid = YamlLoader::loadFromYaml("/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/map_tsudanuma.yaml", 50);
    std::cout << "Loaded grid: " << grid.width << "x" << grid.height << ", resolution: " << grid.resolution << std::endl;
    
    // Generate GVD
    std::cout << "Generating GVD..." << std::endl;
    GvdGenerator gvd;
    auto result = gvd.run(grid);
    
    size_t gvd_count = 0;
    for (auto v : result.gvd_mask) if (v) ++gvd_count;
    std::cout << "GVD pixels: " << gvd_count << std::endl;
    
    // Extract topology with optimized parameters
    std::cout << "Extracting topology with optimized parameters..." << std::endl;
    TopologyExtractor topo;
    
    // Use more aggressive pruning for large maps
    TopologyExtractor::Params params;
    params.prune_min_length = 1.0;  // 1 meter minimum edge length
    params.merge_radius = 0.5;      // 0.5 meter merge radius
    params.resolution = grid.resolution;
    
    topo.setParams(params);
    auto map = topo.run(grid, result.gvd_mask, result.width, result.height, grid.resolution);
    
    std::cout << "Topological map: " << map.nodes.size() << " nodes, " << map.edges.size() << " edges" << std::endl;
    
    // Analyze node distribution
    if (!map.nodes.empty()) {
        double min_x = map.nodes[0].x;
        double max_x = map.nodes[0].x;
        double min_y = map.nodes[0].y;
        double max_y = map.nodes[0].y;
        for (const auto& node : map.nodes) {
            min_x = std::min(min_x, node.x);
            max_x = std::max(max_x, node.x);
            min_y = std::min(min_y, node.y);
            max_y = std::max(max_y, node.y);
        }
        std::cout << "Node bounds: x[" << min_x << ", " << max_x << "], y[" << min_y << ", " << max_y << "]" << std::endl;
        std::cout << "Node range: " << (max_x - min_x) << " x " << (max_y - min_y) << std::endl;
    }
    
    // Save results
    std::cout << "Saving results..." << std::endl;
    
    // Save JSON
    std::ofstream json_file("/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/tsudanuma_new_topo.json");
    json_file << "{\n";
    json_file << "  \"nodes\": [\n";
    for (size_t i = 0; i < map.nodes.size(); ++i) {
        json_file << "    {\"id\": " << map.nodes[i].id << ", \"x\": " << map.nodes[i].x << ", \"y\": " << map.nodes[i].y << "}";
        if (i < map.nodes.size() - 1) json_file << ",";
        json_file << "\n";
    }
    json_file << "  ],\n";
    json_file << "  \"edges\": [\n";
    for (size_t i = 0; i < map.edges.size(); ++i) {
        json_file << "    {\"u\": " << map.edges[i].u << ", \"v\": " << map.edges[i].v << ", \"length\": " << map.edges[i].length << "}";
        if (i < map.edges.size() - 1) json_file << ",";
        json_file << "\n";
    }
    json_file << "  ]\n";
    json_file << "}\n";
    json_file.close();
    
    // Save visualization
    VisualizationOptions options;
    options.image_width = 1200;
    options.image_height = 800;
    options.node_radius = 2;
    options.edge_thickness = 1;
    options.show_node_ids = false;
    options.show_edge_lengths = false;
    
    try {
        Visualizer::saveTopologicalMapAsImage(map, "/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/tsudanuma_new_topo_visual.png", options);
        std::cout << "Visualization saved to tsudanuma_new_topo_visual.png" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Visualization error: " << e.what() << std::endl;
    }
    
    std::cout << "Processing completed successfully!" << std::endl;
    return 0;
}
