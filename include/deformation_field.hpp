#ifndef DEFORMATION_FIELD_HPP
#define DEFORMATION_FIELD_HPP

#include "io.hpp"

#include <eigen3/Eigen/Dense>
#include <cassert>

namespace adi{
    namespace deformation_field{
        class DeformationField{
            public:
                DeformationField(const uint32_t &number_of_grid_points) : m_number_of_grid_points(number_of_grid_points){}

                const double computeScalarPotentialFieldByIndex(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z ); 

                std::vector<Eigen::Vector3d> computeScalarPotentialFieldBasis();
                
            private:
                const uint32_t &m_number_of_grid_points;
        
        };
    }
}
#endif
