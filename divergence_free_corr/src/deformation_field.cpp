#include "deformation_field.hpp"
#include "io.hpp"
#include "constants.hpp"
#include "numerics.hpp"

#include<cmath>
#include<thread>
#include<pthread.h>
#include<algorithm>

namespace adi{
    namespace deformation_field{

        void* VelocityBasisWorker(void *arg)
        {
           VelThreadData* vel_thread_data = static_cast<VelThreadData*>(arg);
            uint32_t thread_start_id = vel_thread_data->start_index;
            uint32_t thread_end_id = vel_thread_data->end_index;
            const std::vector<BasisIndices>& thread_basis_indices = *(vel_thread_data->basis_indices);
            const Eigen::Vector3d& thread_pt = *(vel_thread_data->point);
            DeformationField* thread_df = vel_thread_data->df;
            const Eigen::VectorXd& thread_coeffs_ak = *(vel_thread_data->coeffs_ak);
           
            Eigen::Vector3d* velocity_field = vel_thread_data->point_vel_field;

            for(uint32_t i = thread_start_id;i < thread_end_id;++i)
            {  
                // Compute partial derivatives of the eigenfunction phi
                double coeff_ak;
                if(thread_coeffs_ak.isZero()){
                    coeff_ak = thread_basis_indices.at(i).eigen_value;
                }
                else{
                    coeff_ak = thread_coeffs_ak[i];
                }
                Eigen::Vector3d velocity_field_per_basis_function = thread_df->computeVelocityBasisFunctions(coeff_ak, thread_basis_indices.at(i), thread_pt);
                velocity_field->x() += velocity_field_per_basis_function.x();
                velocity_field->y() += velocity_field_per_basis_function.y();
                velocity_field->z() += velocity_field_per_basis_function.z();   
            }

            return nullptr;
        }

        const double DeformationField::phi(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z )
        {
            return 0.125 * sin(M_PI * index_1 * x) * sin(M_PI * index_2 * y) * sin(M_PI * index_3 * z);
        }

         const double DeformationField::dphidx(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z )
        {
            return 0.5 * M_PI * index_1 * cos(M_PI * index_1 * x) * sin(M_PI * index_2 * y) * sin(M_PI * index_3 * z);
        }

        const double DeformationField::dphidy(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z )
        {
            return 0.5 * M_PI * index_2 * cos(M_PI * index_1 * x) * sin(M_PI * index_2 * y) * sin(M_PI * index_3 * z);
        }

        const double DeformationField::dphidz(const uint32_t &index_1, const uint32_t &index_2,const uint32_t &index_3, const double &x, const double &y, const double &z )
        {
            return 0.5 * M_PI * index_3 * cos(M_PI * index_1 * x) * sin(M_PI * index_2 * y) * sin(M_PI * index_3 * z);
        }

        Eigen::Vector3d DeformationField::computeVelocityBasisFunctions(const double &coeff_ak, const BasisIndices &base_index, const Eigen::Vector3d &point)
        {
            // Compute partial derivatives of the eigenfunction phi
            double dphi_dx = dphidx(base_index.index_1, base_index.index_2, base_index.index_3, point.x(), point.y(), point.z());
            double dphi_dy = dphidy(base_index.index_1, base_index.index_2, base_index.index_3, point.x(), point.y(), point.z());
            double dphi_dz = dphidz(base_index.index_1, base_index.index_2, base_index.index_3, point.x(), point.y(), point.z());
            
            // Compute the basis functions for the eigen mode (eigen basis)
            Eigen::Vector3d v1(0,dphi_dz,-dphi_dy);
            Eigen::Vector3d v2(-dphi_dz, 0, dphi_dx);
            Eigen::Vector3d v3(dphi_dy, -dphi_dx, 0);

            Eigen::Vector3d velocity_at_kth_mode = coeff_ak * (v1 + v2 + v3);

            return velocity_at_kth_mode;
        }

        Eigen::Vector3d DeformationField::computeVelocityField(const Eigen::Vector3d &point, const std::vector<BasisIndices> &base_indices, const Eigen::VectorXd &coeffs_ak)
        {
            Eigen::Vector3d velocity_field;
            velocity_field.setZero();
            
            unsigned int max_num_threads = std::thread::hardware_concurrency();
            pthread_t threads[max_num_threads];
            int thread_ids[max_num_threads];
            VelThreadData vel_thread_data[max_num_threads];

            uint32_t chunk_size = (base_indices.size() + max_num_threads - 1) / (max_num_threads);

            for(size_t t_id = 0;t_id < max_num_threads;++t_id)
            {
                vel_thread_data[t_id].basis_indices = &base_indices;
                vel_thread_data[t_id].coeffs_ak  = &coeffs_ak;
                vel_thread_data[t_id].df = this;
                vel_thread_data[t_id].point = &point;
                vel_thread_data[t_id].start_index = t_id * chunk_size;
                vel_thread_data[t_id].end_index = std::min(static_cast<uint32_t>((t_id + 1) * chunk_size), (uint32_t)base_indices.size());
                vel_thread_data[t_id].point_vel_field = &velocity_field;
                pthread_create(&threads[t_id], NULL,VelocityBasisWorker, &vel_thread_data[t_id]);
            }

            // Join threads
            for (int i = 0; i < max_num_threads; ++i) {
                pthread_join(threads[i], nullptr);
            }
            return velocity_field;
        }

    }
}