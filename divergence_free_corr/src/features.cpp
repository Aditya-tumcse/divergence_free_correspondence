#include "features.hpp"

namespace adi{
    namespace features{
        pcl::PointCloud<pcl::SHOT352>::Ptr computeShotFeatures(const pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals, const double search_radius)
        {
            // Set up SHOT feature estimator
            pcl::SHOTEstimationOMP<pcl::PointNormal, pcl::PointNormal, pcl::SHOT352> shot;
            shot.setInputCloud(cloud_with_normals);
            shot.setInputNormals(cloud_with_normals);
            
            // Create a KdTree for searching
            pcl::search::KdTree<pcl::PointNormal>::Ptr tree(new pcl::search::KdTree<pcl::PointNormal>());
            shot.setSearchMethod(tree);

            // Set the radius for SHOT estimation
            shot.setRadiusSearch(search_radius);

            // Compute SHOT descriptors
            pcl::PointCloud<pcl::SHOT352>::Ptr shot_descriptor(new pcl::PointCloud<pcl::SHOT352>);
            shot.compute(*shot_descriptor);

            return shot_descriptor;
        }
    }
}