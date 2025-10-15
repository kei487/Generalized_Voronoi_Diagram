#pragma once

#include "gvd_topo/core/TopologyExtractor.hpp"
#include <string>

namespace gvd_topo {

struct VisualizationOptions {
    int image_width {800};
    int image_height {600};
    int node_radius {1};
    int edge_thickness {1};
    bool show_node_ids {false};
    bool show_edge_lengths {false};
    double scale_factor {1.0}; // Auto-scale if <= 0
    int background_color[3] {255, 255, 255}; // RGB
    int node_color[3] {0, 255, 0}; // Green
    int edge_color[3] {255, 0, 0}; // Red
    int text_color[3] {0, 0, 0}; // Black
};

class Visualizer {
public:
    static void saveTopologicalMapAsImage(const TopologicalMap& map, 
                                        const std::string& output_path,
                                        const VisualizationOptions& options = VisualizationOptions{});
    
    static void saveTopologicalMapAsImage(const TopologicalMap& map,
                                        const std::string& output_path,
                                        int image_width, int image_height,
                                        double scale_factor = 0.0);
};

} // namespace gvd_topo


