#pragma once

#include <string>

namespace gvd_topo {

struct ProcessingParameters {
    // OccupancyGrid parameters
    int occupancy_threshold {50};
    int morphology_kernel_size {0}; // 0 = disabled
    bool morphology_enabled {false};
    
    // GVD generation parameters
    float distance_epsilon {1e-6f};
    bool use_opencv {true};
    
    // Topology extraction parameters
    double prune_min_length {30.0};
    double merge_radius {30.0};
    int max_trace_steps {100000};
    
    // Output parameters
    std::string output_frame_id {"map"};
    bool publish_debug_info {false};
    
    // Performance parameters
    bool use_parallel_processing {true};
    int max_threads {0}; // 0 = auto-detect
};

struct NodeParameters {
    // ROS topic names
    std::string map_topic {"/map"};
    std::string topological_map_topic {"/topological_map"};
    std::string markers_topic {"/topological_markers"};
    
    // ROS service names
    std::string generate_service {"/generate_topological_map"};
    
    // ROS parameter names (for dynamic_reconfigure)
    std::string param_namespace {"gvd_topo"};
    
    // Processing parameters
    ProcessingParameters processing {};
};

// Helper functions for parameter validation and defaults
ProcessingParameters getDefaultProcessingParameters();
NodeParameters getDefaultNodeParameters();

bool validateParameters(const ProcessingParameters& params);
bool validateParameters(const NodeParameters& params);

} // namespace gvd_topo
