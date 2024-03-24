#include "io.hpp"

int main()
{
    const std::string input_cloud_path = "/home/aditya/Documents/example.ply";
    auto input_cloud = adi::pointCloud(input_cloud_path);
    if(input_cloud.getPointCloudPath().empty())
        return 1;
    else
        return 0;
}