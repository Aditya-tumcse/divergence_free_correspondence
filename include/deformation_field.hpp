#ifndef DEFORMATION_FIELD_HPP
#define DEFORMATION_FIELD_HPP

#include "io.hpp"

#include <eigen3/Eigen/Dense>
#include <cassert>

namespace adi{
    namespace deformation_field{

        struct Indices{
            unsigned int s_x;
            unsigned int s_y;
            unsigned int s_z;

            Indices() : s_x(0), s_y(0), s_z(0){}

            Indices(unsigned int x, unsigned int y, unsigned int z) : s_x(x), s_y(y), s_z(z){}
        };

        const std::vector<std::tuple<unsigned int>> getBasisFunctions(const int nunber_of_basis_functions);
        class deformationField{
            public:
                deformationField(const std::vector<adi::Point> &source_point_cloud, const std::vector<adi::Point> &target_point_cloud, const std::vector<std::pair<adi::Point, adi::Point>> &correspondences);

                Eigen::MatrixXd computeSoftCorrespondences();

                const std::vector<adi::Point> getSourcePointCloud() const { return m_source_point_cloud;}

                const std::vector<adi::Point> getTargetPointCloud() const {return m_target_point_cloud;}

            private:
                const std::vector<adi::Point> &m_source_point_cloud;
                const std::vector<adi::Point> &m_target_point_cloud;
                const std::vector<std::pair<adi::Point, adi::Point>> &m_correspondences;

                const double computeMeanEucledianDistance();
                const double computeMeanDescriptorDistance();
                const Eigen::MatrixXd computeMetricDistance();

                const std::vector<double> computGradientOfBasisFunctions(const std::vector<adi::Point> &point_cloud, const std::vector<Indices> &indices_for_basis_function_computation,const unsigned int &dimension);
        };
    }
}
#endif
