#pragma once

#include <vector>
#include <utility>
#include <cstdint>
#include <string>

namespace gvd_topo {

class OccupancyGrid; // Forward declaration

struct TopoNode {
    int id {0};
    double x {0.0};
    double y {0.0};
    //std::vector<int> edge_ids;
};

struct TopoEdge {
    int id {0};
    int u {0};
    int v {0};
    double length {0.0};
    //std::vector<std::pair<double,double>> polyline; // optional geometry
};

struct TopologicalMap {
    std::vector<TopoNode> nodes;
    std::vector<TopoEdge> edges;
};

class TopologyExtractor {
public:
    struct Params {
        double prune_min_length {0.1};
        double merge_radius {0.05};
        double resolution {0.05};
    };

    TopologyExtractor();
    explicit TopologyExtractor(const Params& p);

    void setParams(const Params& p) { params_ = p; }
    TopologicalMap run(const OccupancyGrid& grid, const std::vector<uint8_t>& gvd_mask, int width, int height, double resolution) const;

private:
    Params params_;
};

// Minimal JSON serialization
std::string toJson(const TopologicalMap& map);

} // namespace gvd_topo


