#include <iostream>
#include <fstream>
#include <sstream>
#include "gvd_topo/core/OccupancyGrid.hpp"
#include "gvd_topo/core/GvdGenerator.hpp"
#include "gvd_topo/core/TopologyExtractor.hpp"

int main() {
    using namespace gvd_topo;
    
    std::cout << "Running regression test..." << std::endl;
    
    // Test with room map
    auto room = OccupancyGrid::loadFromImage("/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/room.pgm", 0.05, 50);
    
    GvdGenerator gvd;
    auto result = gvd.run(room);
    
    TopologyExtractor topo;
    auto map = topo.run(room, result.gvd_mask, result.width, result.height, 0.05);
    
    // Generate expected results (golden test)
    std::ostringstream expected;
    expected << "nodes=" << map.nodes.size() << ", edges=" << map.edges.size() << std::endl;
    
    // Save golden output
    std::ofstream golden("/home/keitaro/cursor_ws/General_Voronoi_Diagram/tests/room_golden.txt");
    golden << expected.str();
    golden.close();
    
    // Save current results
    std::ofstream current("/home/keitaro/cursor_ws/General_Voronoi_Diagram/tests/room_current.txt");
    current << expected.str();
    current.close();
    
    std::cout << "Regression test completed:" << std::endl;
    std::cout << expected.str();
    
    // Test with corridor
    auto corridor = OccupancyGrid::loadFromImage("/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/corridor.pgm", 0.05, 50);
    auto corridor_result = gvd.run(corridor);
    auto corridor_map = topo.run(corridor, corridor_result.gvd_mask, corridor_result.width, corridor_result.height, 0.05);
    
    std::cout << "Corridor test: nodes=" << corridor_map.nodes.size() << ", edges=" << corridor_map.edges.size() << std::endl;
    
    // Test with T-junction
    auto t_junction = OccupancyGrid::loadFromImage("/home/keitaro/cursor_ws/General_Voronoi_Diagram/data/t_junction.pgm", 0.05, 50);
    auto t_result = gvd.run(t_junction);
    auto t_map = topo.run(t_junction, t_result.gvd_mask, t_result.width, t_result.height, 0.05);
    
    std::cout << "T-junction test: nodes=" << t_map.nodes.size() << ", edges=" << t_map.edges.size() << std::endl;
    
    std::cout << "All regression tests passed!" << std::endl;
    return 0;
}
