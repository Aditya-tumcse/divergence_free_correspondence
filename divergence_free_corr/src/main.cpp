#include "../external/include/clue.hpp"
#include "run.hpp"

int main() {
  const std::string source_input_cloud_path =
      "/workspaces/divergence_free_correspondence/data/source_cloud.ply";
  const std::string target_input_cloud_path =
      "/workspaces/divergence_free_correspondence/data/target_cloud.ply";

  LOG_INFO("Starting deformation modelling");
  adi::pointCloud source_cloud(source_input_cloud_path);
  adi::pointCloud target_cloud(target_input_cloud_path);

  run(&source_cloud, &target_cloud);

  return EXIT_SUCCESS;
}