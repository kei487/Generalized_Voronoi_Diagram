#include "gvd_topo/core/GvdGenerator.hpp"
#include "gvd_topo/core/OccupancyGrid.hpp"
#ifdef GVD_TOPO_WITH_OPENCV
#include <opencv2/imgproc.hpp>
#include <opencv2/core.hpp>
#endif

namespace gvd_topo {

GvdGenerator::GvdGenerator() = default;
GvdGenerator::GvdGenerator(const Params& p) : params_(p) {}

GvdResult GvdGenerator::run(const OccupancyGrid& grid) const {
    GvdResult result;
    result.width = grid.width;
    result.height = grid.height;
    result.distance.assign(static_cast<size_t>(grid.width * grid.height), 0.0f);
    result.gvd_mask.assign(static_cast<size_t>(grid.width * grid.height), 0);
    // Simple EDT via OpenCV if available
#ifdef GVD_TOPO_WITH_OPENCV
    if (!grid.empty()) {
        cv::Mat occ(grid.height, grid.width, CV_8UC1);
        for (int y = 0; y < grid.height; ++y) {
            uint8_t* row = occ.ptr<uint8_t>(y);
            for (int x = 0; x < grid.width; ++x) {
                int8_t c = grid.data[grid.index(x, y)];
                row[x] = (c == static_cast<int8_t>(100)) ? 0 : 255;
            }
        }
        cv::Mat dist;
        cv::distanceTransform(occ, dist, cv::DIST_L2, 3);
        // copy back
        for (int y = 0; y < grid.height; ++y) {
            const float* row = dist.ptr<float>(y);
            for (int x = 0; x < grid.width; ++x) {
                result.distance[grid.index(x, y)] = row[x] * static_cast<float>(grid.resolution);
            }
        }

        // Morphological skeletonization (similar to scikit-image's medial_axis)
        // Step 1: Threshold distance map to create binary image (free space)
        const int w = grid.width;
        const int h = grid.height;
        cv::Mat binary(h, w, CV_8UC1);
        const float dist_threshold = 0.001f; // cells with distance > threshold are considered free
        
        for (int y = 0; y < h; ++y) {
            uint8_t* row = binary.ptr<uint8_t>(y);
            for (int x = 0; x < w; ++x) {
                row[x] = (result.distance[y * w + x] > dist_threshold) ? 255 : 0;
            }
        }
        
        // Step 2: Morphological skeletonization using iterative thinning
        // This approximates scikit-image's medial_axis behavior
        cv::Mat skeleton = cv::Mat::zeros(h, w, CV_8UC1);
        cv::Mat temp;
        cv::Mat eroded;
        
        cv::Mat element = cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3));
        
        bool done = false;
        int iteration = 0;
        const int max_iterations = 3000;//std::min(std::max(w, h), 10000); // Limit based on map size
        while (!done && iteration < max_iterations) {
            cv::erode(binary, eroded, element);
            cv::dilate(eroded, temp, element);
            cv::subtract(binary, temp, temp);
            cv::bitwise_or(skeleton, temp, skeleton);
            eroded.copyTo(binary);
            
            done = (cv::countNonZero(binary) == 0);
            ++iteration;
        }


        // Step 3: Copy skeleton to GVD mask
        for (int y = 0; y < h; ++y) {
            const uint8_t* row = skeleton.ptr<uint8_t>(y);
            for (int x = 0; x < w; ++x) {
                result.gvd_mask[y * w + x] = row[x];
            }
        }
    }
#endif
    return result;
}

} // namespace gvd_topo


