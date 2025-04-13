#include "../external/include/clue.hpp"
#include "run.hpp"
#include <filesystem>

int main(int argc, char **argv) {
  std::string data_folder = " ";
  if (argc > 1) {
    data_folder = argv[1];
    if (!data_folder.empty() && data_folder.back() != '/')
      data_folder += '/';
  }
  const std::string source_input_cloud_path = data_folder + "source_cloud.ply";
  const std::string target_input_cloud_path = data_folder + "target_cloud.ply";

  if (std::filesystem::exists(source_input_cloud_path) &&
      std::filesystem::exists(target_input_cloud_path)) {
    LOG_INFO("Starting deformation modelling");
    adi::pointCloud source_cloud(source_input_cloud_path);
    adi::pointCloud target_cloud(target_input_cloud_path);

    run(&source_cloud, &target_cloud);

    return EXIT_SUCCESS;
  } else {
    LOG_INFO("No input files found. Hence exiting!");
    return EXIT_FAILURE;
  }
}