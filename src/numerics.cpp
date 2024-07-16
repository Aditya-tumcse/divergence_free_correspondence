#include <random>

#include "numerics.hpp"
#include "deformation_field.hpp"

namespace numerics{
    const std::vector<adi::deformation_field::BasisIndices> GenerateBasisIndices(const uint32_t &max_number_of_velocity_basis)
    {
        std::vector<adi::deformation_field::BasisIndices> basis_indices;
        for(uint32_t i = 1;i <= max_number_of_velocity_basis;++i)
        {
            for(uint32_t j = 1;j <= max_number_of_velocity_basis;++j)
            {
                for(uint32_t k = 1;k <= max_number_of_velocity_basis;++k)
                {
                    adi::deformation_field::BasisIndices base_index;
                    base_index.index_1 = i;
                    base_index.index_2 = j;
                    base_index.index_3 = k;
                    double eigen_val = -std::pow(M_PI,2) * (std::pow(i,2) + std::pow(j,2) + std::pow(k,2));
                    base_index.eigen_value = eigen_val;
                    
                    basis_indices.emplace_back(base_index);
                }
            }
        }

        return basis_indices;
    }

    const double ObtainCoefficientOfVelocityField(const double eigen_value)
    {
        const double mean = 0.0;
        const double stddev = eigen_value;

        // Create a random number generator and a normal distribution object
        std::random_device rd;  // Random device for seeding
        std::mt19937 gen(rd()); // Mersenne Twister generator seeded with rd()
        std::normal_distribution<> d(mean, stddev);

        // Generate a sample from the normal distribution
        double sample = d(gen);

        return sample;
    }

    Eigen::Vector3d RungeKutaIntegration(const Eigen::Vector3d &src_pt, const double dt)
    {
        auto basis_indices = GenerateBasisIndices(MAX_NUMBER_OF_VELOCITY_BASIS);
        auto df = adi::deformation_field::DeformationField();

        Eigen::Vector3d k1 = df.computeVelocityField(src_pt, basis_indices);
        Eigen::Vector3d midpoint = src_pt + 0.5 * dt * k1;
        Eigen::Vector3d k2 = df.computeVelocityField(midpoint, basis_indices);

        return src_pt + dt * k2;
    }
}