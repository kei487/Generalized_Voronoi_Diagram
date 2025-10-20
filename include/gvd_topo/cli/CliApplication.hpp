#pragma once

#include "gvd_topo/utils/ConfigManager.hpp"
#include "gvd_topo/core/OccupancyGrid.hpp"
#include "gvd_topo/core/GvdGenerator.hpp"
#include "gvd_topo/core/TopologyExtractor.hpp"
#include "gvd_topo/core/Visualizer.hpp"
#include "gvd_topo/utils/Timer.hpp"
#include <memory>

namespace gvd_topo {

/**
 * @brief Command Line Interface Application
 * 
 * This class encapsulates the main application logic for the CLI tool,
 * separating it from the command line parsing and configuration management.
 */
class CliApplication {
public:
    CliApplication();
    ~CliApplication() = default;

    /**
     * @brief Run the application with given configuration
     * @param config Configuration to use
     * @return Exit code (0 for success, non-zero for error)
     */
    int run(const ConfigOptions& config);

    /**
     * @brief Set timing report callback
     * @param callback Function to call for timing reports
     */
    void setTimingCallback(std::function<void(const std::string&, double)> callback);

private:
    std::unique_ptr<GvdGenerator> gvd_generator_;
    std::unique_ptr<TopologyExtractor> topology_extractor_;
    std::function<void(const std::string&, double)> timing_callback_;

    /**
     * @brief Load occupancy grid based on configuration
     * @param config Configuration containing input parameters
     * @return Loaded occupancy grid
     */
    OccupancyGrid loadOccupancyGrid(const ConfigOptions& config);

    /**
     * @brief Generate GVD from occupancy grid
     * @param grid Input occupancy grid
     * @return GVD result
     */
    GvdResult generateGVD(const OccupancyGrid& grid);

    /**
     * @brief Extract topology from GVD result
     * @param grid Input occupancy grid
     * @param gvd_result GVD result
     * @param resolution Map resolution
     * @return Topological map
     */
    TopologicalMap extractTopology(const OccupancyGrid& grid, const GvdResult& gvd_result, double resolution);

    /**
     * @brief Save outputs based on configuration
     * @param config Configuration containing output parameters
     * @param grid Input occupancy grid
     * @param gvd_result GVD result
     * @param topo_map Topological map
     */
    void saveOutputs(const ConfigOptions& config, 
                    const OccupancyGrid& grid,
                    const GvdResult& gvd_result,
                    const TopologicalMap& topo_map);

    /**
     * @brief Print processing statistics
     * @param gvd_result GVD result
     * @param topo_map Topological map
     */
    void printStatistics(const GvdResult& gvd_result, const TopologicalMap& topo_map);
};

} // namespace gvd_topo

