#include <chrono>

#include "run.hpp"

void run(adi::pointCloud *source_cloud, adi::pointCloud *target_cloud)
{
    // downsample the point clouds
    source_cloud->samplePointCloud(NUMBER_OF_SAMPLE_POINTS);
    target_cloud->samplePointCloud(NUMBER_OF_SAMPLE_POINTS);
   
    // Compute features
    source_cloud->computeShotFeatures(SEARCH_RADIUS);
    std::cout << "Completed computing SHOT features for source cloud " << std::endl;

    target_cloud->computeShotFeatures(SEARCH_RADIUS);
    std::cout << "Completed computing SHOT features for target cloud " << std::endl;

    // get downsampled cloud with features
    std::vector<adi::Point> source_downsampled_cloud = source_cloud->getPointCloud();
    std::vector<adi::Point> target_downsampled_cloud = target_cloud->getPointCloud();

    std::cout << "Size of source cloud for processing: " << source_downsampled_cloud.size() << std::endl;
    std::cout << "Size of target clou for processing: " << target_downsampled_cloud.size() << std::endl;

    // Compute initial correspondence based on SHOT features
    std::unique_ptr<std::vector<std::pair<adi::Point, adi::Point>>> initial_correspondences = utilities::computeCorrespondences(source_downsampled_cloud, target_downsampled_cloud);

    // Gauss Newton optimization looop until convergence
    uint32_t iter = 0;
    std::vector<adi::deformation_field::BasisIndices> base_indices = adi::numerics::GenerateBasisIndices(MAX_NUMBER_OF_VELOCITY_BASIS);

    //TODO:precompute the velocity basis functions as a matrix


    Eigen::VectorXd coeffs_ak;
    coeffs_ak.resize(base_indices.size());
    coeffs_ak.setZero();
    while(iter < MAX_NUMBER_OF_ITERS)
    {
        // Update point cloud Y with the current deformation field
        std::vector<adi::Point> updated_cloud;
        updated_cloud.reserve(source_downsampled_cloud.size());

        std::cout << "Starting RK2 integration" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        // for(uint32_t i = 0; i < source_downsampled_cloud.size(); ++i) {
        //     Eigen::Vector3d updated_pt =  adi::numerics::RungeKutaIntegration(source_downsampled_cloud.at(i).s_point, base_indices,coeffs_ak, NUMBER_OF_TIME_STEPS); // Adjust dt as needed
        //     updated_cloud[i].s_point.x() = updated_pt.x();
        //     updated_cloud[i].s_point.y() = updated_pt.y();
        //     updated_cloud[i].s_point.z() = updated_pt.z();        
        // }

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
        std::cout << "Time taken for RK integration for point ID: " << duration.count() << std::endl;
        std::cout << "Completed RK integration" << std::endl;

        // E-step : Compute soft correspondences
        auto matches = adi::matching::Matching(*initial_correspondences);
        auto soft_corr_matrix = matches.computeSoftCorrespondences(updated_cloud, target_downsampled_cloud);

        // M-step : Gauss Newton optimization
        Eigen::MatrixXd hessian_matrix = adi::numerics::ComputeHessian(base_indices, updated_cloud, target_downsampled_cloud, soft_corr_matrix, r0);

        Eigen::VectorXd gradient = adi::numerics::ComputeGradient(base_indices, updated_cloud, target_downsampled_cloud, soft_corr_matrix, r0);

        Eigen::VectorXd delta_a = -hessian_matrix.ldlt().solve(gradient);

        coeffs_ak += delta_a;

        std::cout << "Updated coefficients of deformation field for iteration " << iter << std::endl;
        if(delta_a.norm() < TOLERANCE)
        {
            std::cout << "Coefficients of deformation field converged" << std::endl;
            break;
        }
        iter = iter + 1;

    }

    for(size_t a_k_i = 0;a_k_i < coeffs_ak.size();++a_k_i)
    {
        std::cout << coeffs_ak[a_k_i] << std::endl;
    }
    return;   
}