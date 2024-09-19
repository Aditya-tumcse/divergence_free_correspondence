#include "io.hpp"
#include "constants.hpp"
#include "deformation_field.hpp"
#include "numerics.hpp"
#include "matching.hpp"

int main()
{
    const std::string source_input_cloud_path = "/home/aditya/Downloads/MPI-FAUST/training/registrations/tr_reg_000.ply";
    const std::string target_input_cloud_path = "/home/aditya/Downloads/MPI-FAUST/training/registrations/tr_reg_002.ply";
    
    // get the source and target point cloud
    auto source_input_cloud = adi::pointCloud(source_input_cloud_path, SEARCH_RADIUS);
    std::cout << "Size of input point cloud: " << source_input_cloud.getSize() << std::endl;
    auto target_input_cloud = adi::pointCloud(target_input_cloud_path, SEARCH_RADIUS);
    std::cout << "Size of target point cloud: " << target_input_cloud.getSize() << std::endl;

    // downsample the point clouds
    auto source_downsampled_cloud = source_input_cloud.samplePointCloud(NUMBER_OF_SAMPLE_POINTS);
    auto target_downsampled_cloud = target_input_cloud.samplePointCloud(NUMBER_OF_SAMPLE_POINTS);
    
    // Compute initial correspondence based on SHOT features
    std::unique_ptr<std::vector<std::pair<adi::Point, adi::Point>>> initial_correspondences = utilities::computeCorrespondences(source_input_cloud.getPointCloud(), target_input_cloud.getPointCloud());

    // Gauss Newton optimization looop until convergence
    uint32_t iter = 0;
    std::vector<adi::deformation_field::BasisIndices> base_indices = adi::numerics::GenerateBasisIndices(MAX_NUMBER_OF_VELOCITY_BASIS);
    Eigen::VectorXd coeffs_ak;
    coeffs_ak.resize(base_indices.size());
    coeffs_ak.setZero();
    while(iter < MAX_NUMBER_OF_ITERS)
    {
        // Update point cloud Y with the current deformation field
        std::vector<adi::Point> updated_cloud;
        updated_cloud.reserve(source_downsampled_cloud->size());
        for(uint32_t i = 0; i < source_downsampled_cloud->size(); ++i) {
            Eigen::Vector3d updated_pt =  adi::numerics::RungeKutaIntegration(source_downsampled_cloud->at(i).s_point, base_indices,coeffs_ak, 1000); // Adjust dt as needed
            updated_cloud[i].s_point.x() = updated_pt.x();
            updated_cloud[i].s_point.y() = updated_pt.y();
            updated_cloud[i].s_point.z() = updated_pt.z();        
        }
        
        // E-step : Compute soft correspondences
        auto matches = adi::matching::Matching(*initial_correspondences);
        auto soft_corr_matrix = matches.computeSoftCorrespondences(updated_cloud, *target_downsampled_cloud);

        // M-step : Gauss Newton optimization
        Eigen::MatrixXd hessian_matrix = adi::numerics::ComputeHessian(base_indices, updated_cloud, *target_downsampled_cloud, soft_corr_matrix, r0);

        Eigen::VectorXd gradient = adi::numerics::ComputeGradient(base_indices, updated_cloud, *target_downsampled_cloud, soft_corr_matrix, r0);

        Eigen::VectorXd delta_a = -hessian_matrix.ldlt().solve(gradient);

        coeffs_ak += delta_a;

        if(delta_a.norm() < TOLERANCE)
            break;

    }
    return 0;

}