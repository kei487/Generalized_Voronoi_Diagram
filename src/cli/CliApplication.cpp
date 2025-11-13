#include "gvd_topo/cli/CliApplication.hpp"
#include "gvd_topo/io/YamlLoader.hpp"
#include <iostream>
#include <fstream>

#ifdef GVD_TOPO_WITH_OPENCV
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#endif

namespace gvd_topo {

CliApplication::CliApplication() 
    : gvd_generator_(std::make_unique<GvdGenerator>())
    , topology_extractor_(std::make_unique<TopologyExtractor>())
    , timing_callback_([](const std::string&, double){}) {
}

void CliApplication::setTimingCallback(std::function<void(const std::string&, double)> callback) {
    timing_callback_ = callback;
}

int CliApplication::run(const ConfigOptions& config) {
    try {
        // Load occupancy grid
        OccupancyGrid grid;
        {
            ScopeTimer timer("load+preprocess", timing_callback_);
            grid = loadOccupancyGrid(config);
        }

        // Generate GVD
        GvdResult gvd_result;
        {
            ScopeTimer timer("EDT+GVD", timing_callback_);
            gvd_result = generateGVD(grid);
        }

        // Extract topology
        TopologicalMap topo_map;
        {
            ScopeTimer timer("topology", timing_callback_);
            topo_map = extractTopology(grid, gvd_result, config.input.resolution);
        }

        // Print statistics
        printStatistics(gvd_result, topo_map);

        // Save outputs
        saveOutputs(config, grid, gvd_result, topo_map);

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

OccupancyGrid CliApplication::loadOccupancyGrid(const ConfigOptions& config) {
    if (config.benchmark.enabled) {
        return OccupancyGrid::randomMap(
            config.benchmark.width, 
            config.benchmark.height, 
            config.input.resolution, 
            config.benchmark.occupancy_ratio, 
            config.benchmark.seed
        );
    } else if (!config.input.yaml_file.empty()) {
        return YamlLoader::loadFromYaml(config.input.yaml_file, config.input.occupancy_threshold);
    } else if (!config.input.image_file.empty()) {
        return OccupancyGrid::loadFromImage(config.input.image_file, config.input.resolution, config.input.occupancy_threshold);
    } else {
        // Default empty grid
        return OccupancyGrid(10, 10, config.input.resolution);
    }
}

GvdResult CliApplication::generateGVD(const OccupancyGrid& grid) {
    return gvd_generator_->run(grid);
}

TopologicalMap CliApplication::extractTopology(const OccupancyGrid& grid, const GvdResult& gvd_result, double resolution) {
    return topology_extractor_->run(grid, gvd_result.gvd_mask, gvd_result.width, gvd_result.height, resolution);
}

void CliApplication::saveOutputs(const ConfigOptions& config, 
                                const OccupancyGrid& grid,
                                const GvdResult& gvd_result,
                                const TopologicalMap& topo_map) {
    // Save topological map JSON
    if (!config.output.map_file.empty()) {
        std::ofstream ofs(config.output.map_file);
        if (!ofs) {
            std::cerr << "Failed to open " << config.output.map_file << std::endl;
            return;
        }
        ofs << toJson(topo_map) << std::endl;
        std::cout << "Wrote map: " << config.output.map_file << std::endl;
    }

#ifdef GVD_TOPO_WITH_OPENCV
    // Save GVD visualization
    if (!config.output.gvd_image.empty()) {
        ScopeTimer timer("overlay", timing_callback_);
        cv::Mat vis(grid.height, grid.width, CV_8UC3);
        
        // Create base visualization
        for (int y = 0; y < grid.height; ++y) {
            for (int x = 0; x < grid.width; ++x) {
                int8_t c = grid.data[grid.index(x, y)];
                uint8_t v = 127;
                if (c == static_cast<int8_t>(100)) v = 0; // occupied -> black
                else if (c == static_cast<int8_t>(0)) v = 255; // free -> white
                vis.at<cv::Vec3b>(y,x) = cv::Vec3b(v,v,v);
            }
        }
        
        // Overlay GVD
        std::cout << "overlay gvd" << std::endl;
        for (int y = 0; y < gvd_result.height; ++y) {
            for (int x = 0; x < gvd_result.width; ++x) {
                if (gvd_result.gvd_mask[y * gvd_result.width + x]) {
                    vis.at<cv::Vec3b>(y,x) = cv::Vec3b(255,0,0); // GVD -> red
                }
            }
        }
        
        // Draw edges
        // for (const auto& edge : topo_map.edges) {
        //     for (size_t i = 1; i < edge.polyline.size(); ++i) {
        //         cv::Point p0(static_cast<int>(std::round(edge.polyline[i-1].first / config.input.resolution)),
        //                      static_cast<int>(std::round(edge.polyline[i-1].second / config.input.resolution)));
        //         cv::Point p1(static_cast<int>(std::round(edge.polyline[i].first / config.input.resolution)),
        //                      static_cast<int>(std::round(edge.polyline[i].second / config.input.resolution)));
        //         cv::line(vis, p0, p1, cv::Scalar(0,0,255), 1); // edges -> blue
        //     }
        // }
        
        // Draw nodes
        // for (const auto& node : topo_map.nodes) {
        //     cv::Point p(static_cast<int>(std::round(node.x / config.input.resolution)),
        //                 static_cast<int>(std::round(node.y / config.input.resolution)));
        //     cv::circle(vis, p, 2, cv::Scalar(0,255,0), -1); // nodes -> green
        // }
        
        if (!cv::imwrite(config.output.gvd_image, vis)) {
            std::cerr << "Failed to write " << config.output.gvd_image << std::endl;
            return;
        }
        std::cout << "Wrote overlay: " << config.output.gvd_image << std::endl;
    }
#endif

    // Save topological map visualization
    if (!config.output.topo_image.empty()) {
        ScopeTimer timer("topo-visualization", timing_callback_);
        try {
            Visualizer::saveTopologicalMapAsImage(topo_map, config.output.topo_image, gvd_result.width, gvd_result.height);
            std::cout << "Wrote topological map image: " << config.output.topo_image << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Failed to create topological map image: " << e.what() << std::endl;
        }
    }
}

void CliApplication::printStatistics(const GvdResult& gvd_result, const TopologicalMap& topo_map) {
    size_t gvd_count = 0;
    for (auto v : gvd_result.gvd_mask) {
        if (v) ++gvd_count;
    }

    std::cout << "gvd_topo_cli: (" << gvd_result.width << "x" << gvd_result.height << ")" << std::endl;
    std::cout << "GVD pixels=" << gvd_count << std::endl;
    std::cout << "nodes=" << topo_map.nodes.size() << ", edges=" << topo_map.edges.size() << std::endl;
}

} // namespace gvd_topo

