#include "gvd_topo/core/Visualizer.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#ifdef GVD_TOPO_WITH_OPENCV
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/core.hpp>
#endif

namespace gvd_topo {

void Visualizer::saveTopologicalMapAsImage(const TopologicalMap& map, 
                                          const std::string& output_path,
                                          const VisualizationOptions& options) {
#ifdef GVD_TOPO_WITH_OPENCV
    if (map.nodes.empty()) {
        std::cerr << "Warning: No nodes to visualize" << std::endl;
        return;
    }
    
    // Calculate bounds
    double min_x = map.nodes[0].x, max_x = map.nodes[0].x;
    double min_y = map.nodes[0].y, max_y = map.nodes[0].y;
    
    for (const auto& node : map.nodes) {
        min_x = std::min(min_x, node.x);
        max_x = std::max(max_x, node.x);
        min_y = std::min(min_y, node.y);
        max_y = std::max(max_y, node.y);
    }
    
    // Add margin (at least 10% or 1.0 unit)
    double range_x = max_x - min_x;
    double range_y = max_y - min_y;
    double margin_x = std::max(range_x * 0.1, 1.0);
    double margin_y = std::max(range_y * 0.1, 1.0);
    min_x -= margin_x;
    max_x += margin_x;
    min_y -= margin_y;
    max_y += margin_y;
    
    // Calculate scale
    double scale;
    if (options.scale_factor > 0) {
        scale = options.scale_factor;
    } else {
        double scale_x = static_cast<double>(options.image_width) / (max_x - min_x);
        double scale_y = static_cast<double>(options.image_height) / (max_y - min_y);
        scale = std::min(scale_x, scale_y);
    }
    
    // Create image
    cv::Mat image(options.image_height, options.image_width, CV_8UC3);
    image.setTo(cv::Scalar(options.background_color[0], options.background_color[1], options.background_color[2]));
    
    // Draw edges first
    for (const auto& edge : map.edges) {
        if (edge.polyline.empty()) continue;
        
        std::vector<cv::Point> points;
        for (const auto& pt : edge.polyline) {
            int x = static_cast<int>((pt.first - min_x) * scale);
            int y = static_cast<int>((pt.second - min_y) * scale);
            y = options.image_height - y; // Flip Y coordinate
            points.push_back(cv::Point(x, y));
        }
        
        if (points.size() > 1) {
            for (size_t i = 1; i < points.size(); ++i) {
                cv::line(image, points[i-1], points[i], 
                        cv::Scalar(options.edge_color[0], options.edge_color[1], options.edge_color[2]), 
                        options.edge_thickness);
            }
        }
        
        // Draw edge length if requested
        if (options.show_edge_lengths && !edge.polyline.empty()) {
            const auto& mid_pt = edge.polyline[edge.polyline.size() / 2];
            int x = static_cast<int>((mid_pt.first - min_x) * scale);
            int y = static_cast<int>((mid_pt.second - min_y) * scale);
            y = options.image_height - y;
            
            std::string length_text = std::to_string(static_cast<int>(edge.length * 100) / 100.0);
            cv::putText(image, length_text, cv::Point(x, y), 
                       cv::FONT_HERSHEY_SIMPLEX, 0.3, 
                       cv::Scalar(options.text_color[0], options.text_color[1], options.text_color[2]), 1);
        }
    }
    
    // Draw nodes
    for (const auto& node : map.nodes) {
        int x = static_cast<int>((node.x - min_x) * scale);
        int y = static_cast<int>((node.y - min_y) * scale);
        y = options.image_height - y; // Flip Y coordinate
        
        cv::circle(image, cv::Point(x, y), 0.5, //options.node_radius, 
                  cv::Scalar(options.node_color[0], options.node_color[1], options.node_color[2]), -1);
        
        // Draw node ID if requested
        if (options.show_node_ids) {
            std::string id_text = std::to_string(node.id);
            cv::putText(image, id_text, cv::Point(x + options.node_radius + 2, y), 
                       cv::FONT_HERSHEY_SIMPLEX, 0.4, 
                       cv::Scalar(options.text_color[0], options.text_color[1], options.text_color[2]), 1);
        }
    }
    
    // Save image
    if (!cv::imwrite(output_path, image)) {
        throw std::runtime_error("Failed to save topological map image: " + output_path);
    }
#else
    (void)map; (void)output_path; (void)options;
    throw std::runtime_error("OpenCV is required for topological map visualization");
#endif
}

void Visualizer::saveTopologicalMapAsImage(const TopologicalMap& map,
                                          const std::string& output_path,
                                          int image_width, int image_height,
                                          double scale_factor) {
    VisualizationOptions options;
    options.image_width = image_width;
    options.image_height = image_height;
    options.scale_factor = scale_factor;
    saveTopologicalMapAsImage(map, output_path, options);
}

} // namespace gvd_topo
