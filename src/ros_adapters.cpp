#include "gvd_topo/ros_adapters.hpp"
#include <cmath>

namespace gvd_topo {

OccupancyGrid fromMsg(const OccupancyGridMsg& msg) {
    OccupancyGrid grid;
    grid.width = msg.width;
    grid.height = msg.height;
    grid.resolution = msg.resolution;
    grid.origin.x = msg.origin.x;
    grid.origin.y = msg.origin.y;
    grid.origin.theta = msg.origin.theta;
    grid.data = msg.data;
    return grid;
}

OccupancyGridMsg toMsg(const OccupancyGrid& grid) {
    OccupancyGridMsg msg;
    msg.width = grid.width;
    msg.height = grid.height;
    msg.resolution = grid.resolution;
    msg.origin.x = grid.origin.x;
    msg.origin.y = grid.origin.y;
    msg.origin.theta = grid.origin.theta;
    msg.data = grid.data;
    return msg;
}

MarkerArray topologicalMapToMarkers(const TopologicalMap& map, const std::string& frame_id) {
    MarkerArray array;
    (void)frame_id; // placeholder for future ROS frame_id usage
    
    // Add node markers (spheres)
    for (const auto& node : map.nodes) {
        Marker marker;
        marker.id = node.id;
        marker.pose.x = node.x;
        marker.pose.y = node.y;
        marker.scale = 0.5;
        marker.color = {0, 255, 0, 255}; // Green for nodes
        array.markers.push_back(marker);
    }
    
    // Add edge markers (lines)
    for (const auto& edge : map.edges) {
        Marker marker;
        marker.id = edge.id + 10000; // Offset to avoid ID conflicts
        marker.pose.x = 0.0; //(edge.polyline.empty()) ? 0.0 : edge.polyline[0].first;
        marker.pose.y = 0.0; //(edge.polyline.empty()) ? 0.0 : edge.polyline[0].second;
        marker.scale = edge.length;
        marker.color = {0, 0, 255, 255}; // Blue for edges
        array.markers.push_back(marker);
    }
    
    return array;
}

std::vector<Marker> nodesToMarkers(const std::vector<TopoNode>& nodes, const std::string& frame_id) {
    (void)frame_id;
    std::vector<Marker> markers;
    markers.reserve(nodes.size());
    
    for (const auto& node : nodes) {
        Marker marker;
        marker.id = node.id;
        marker.pose.x = node.x;
        marker.pose.y = node.y;
        marker.scale = 0.3;
        marker.color = {0, 255, 0, 255}; // Green for nodes
        markers.push_back(marker);
    }
    
    return markers;
}

std::vector<Marker> edgesToMarkers(const std::vector<TopoEdge>& edges, const std::string& frame_id) {
    (void)frame_id;
    std::vector<Marker> markers;
    markers.reserve(edges.size());
    
    for (const auto& edge : edges) {
        Marker marker;
        marker.id = edge.id + 10000; // Offset to avoid ID conflicts
        marker.pose.x = 0.0;//(edge.polyline.empty()) ? 0.0 : edge.polyline[0].first;
        marker.pose.y = 0.0;//(edge.polyline.empty()) ? 0.0 : edge.polyline[0].second;
        marker.scale = edge.length;
        marker.color = {255, 0, 0, 255}; // Red for edges
        markers.push_back(marker);
    }
    
    return markers;
}

} // namespace gvd_topo
