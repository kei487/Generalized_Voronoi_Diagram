#include <iostream>
#include "gvd_topo/ros_adapters.hpp"
#include "gvd_topo/parameters.hpp"

int main() {
    using namespace gvd_topo;
    
    // Test ROS adapters
    std::cout << "Testing ROS adapters..." << std::endl;
    
    // Create a test occupancy grid
    OccupancyGrid grid(5, 5, 0.1);
    for (int i = 0; i < 25; ++i) {
        grid.data[i] = (i % 3 == 0) ? static_cast<int8_t>(100) : static_cast<int8_t>(0);
    }
    
    // Convert to ROS message format
    OccupancyGridMsg msg = toMsg(grid);
    std::cout << "OccupancyGridMsg: " << msg.width << "x" << msg.height << " @ " << msg.resolution << std::endl;
    
    // Convert back
    OccupancyGrid grid2 = fromMsg(msg);
    std::cout << "Round-trip conversion: " << (grid2.width == grid.width ? "OK" : "FAIL") << std::endl;
    
    // Test marker generation
    TopologicalMap map;
    TopoNode node1; node1.id = 1; node1.x = 1.0; node1.y = 2.0;
    TopoNode node2; node2.id = 2; node2.x = 3.0; node2.y = 4.0;
    map.nodes = {node1, node2};
    
    TopoEdge edge; edge.id = 1; edge.u = 1; edge.v = 2; edge.length = 2.828;
    //edge.polyline = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}};
    map.edges = {edge};
    
    MarkerArray markers = topologicalMapToMarkers(map);
    std::cout << "Generated " << markers.markers.size() << " markers" << std::endl;
    
    // Test parameters
    std::cout << "Testing parameters..." << std::endl;
    ProcessingParameters proc_params = getDefaultProcessingParameters();
    NodeParameters node_params = getDefaultNodeParameters();
    
    std::cout << "Default occupancy threshold: " << proc_params.occupancy_threshold << std::endl;
    std::cout << "Default map topic: " << node_params.map_topic << std::endl;
    std::cout << "Parameter validation: " << (validateParameters(proc_params) ? "OK" : "FAIL") << std::endl;
    
    std::cout << "ROS adapters test completed successfully!" << std::endl;
    return 0;
}
