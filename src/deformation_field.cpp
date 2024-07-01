#include "deformation_field.hpp"
#include "io.hpp"
#include "constants.hpp"

#include<cmath>

namespace adi{
    namespace deformation_field{
        const double DeformationField::computeScalarPotentialFieldByIndex(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z )
        {
            return 0.5 * sin(M_PI * index_1 * x) * sin(M_PI * index_2 * y) * sin(M_PI * index_3 * z);
        }


    }
}