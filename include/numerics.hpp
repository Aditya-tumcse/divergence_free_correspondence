#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include<eigen3/Eigen/Dense>
#include "deformation_field.hpp"
#include "constants.hpp"

namespace numerics{
    const std::vector<adi::deformation_field::BasisIndices> GenerateBasisIndices(const uint32_t &max_number_of_velocity_basis);

    const double ObtainCoefficientOfVelocityField(const double eigen_value);

    Eigen::Vector3d RungeKutaIntegration(const Eigen::Vector3d &src_pt, const double dt);
}
#endif