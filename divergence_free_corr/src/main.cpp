#include "run.hpp"

int main()
{
    const std::string source_input_cloud_path = "../data/source_cloud.ply";
    const std::string target_input_cloud_path = "../data/target_cloud.ply";
    
    run(source_input_cloud_path, target_input_cloud_path);

    return EXIT_SUCCESS;
}