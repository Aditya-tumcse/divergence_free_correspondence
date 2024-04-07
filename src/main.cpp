#include "io.hpp"

int main()
{
    const std::string source_input_cloud_path = "/home/aditya/Documents/example.ply";
    const std::string target_input_cloud_path = "";
    
    // get the source and target point cloud
    auto source_input_cloud = adi::pointCloud(source_input_cloud_path);
    auto target_input_cloud = adi::pointCloud(target_input_cloud_path);
    
    // if either of the pointcloud paths are empty, exit with failure
    if(source_input_cloud.getPointCloudPath().empty() || target_input_cloud.getPointCloudPath().empty())
        return 1;
    else
        return 0;
}