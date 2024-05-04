#include "io.hpp"
#include "constants.hpp"
#include "deformation_field.hpp"

int main()
{
    const std::string source_input_cloud_path = "/home/aditya/Documents/example.ply";
    const std::string target_input_cloud_path = "";
    
    // get the source and target point cloud
    auto source_input_cloud = adi::pointCloud(source_input_cloud_path, SEARCH_RADIUS);
    auto target_input_cloud = adi::pointCloud(target_input_cloud_path, SEARCH_RADIUS);

    // downsample the point clouds
    source_input_cloud.samplePointCloud(NUMBER_OF_SAMPLE_POINTS);
    target_input_cloud.samplePointCloud(NUMBER_OF_SAMPLE_POINTS);
    
    // Compute initial correspondence based on SHOT features
    std::vector<std::pair<adi::Point, adi::Point>> initial_correspondences = utilities::computeCorrespondences(source_input_cloud.getPointCloud(), target_input_cloud.getPointCloud());

    //Initialize deformation field
    adi::deformation_field::deformationField df(initial_correspondences);

    return 0;

}