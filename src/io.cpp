#include "io.hpp"

namespace adi{
    pcl::PointCloud<pcl::PointXYZ>::Ptr pointCloud::readPointCloud(const std::string &point_cloud_path)
    {
        if(point_cloud_path.empty())
            return nullptr;
        else{
            pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud(new pcl::PointCloud<pcl::PointXYZ>);
            if(pcl::io::loadPLYFile<pcl::PointXYZ>(point_cloud_path, *input_cloud) == -1){
                std::cerr << "Error loading point cloud file" << std::endl;
                return nullptr;
            }
            else{   
                return input_cloud;
            }
        }
    }

    pcl::PointCloud<pcl::PointNormal>::Ptr pointCloud::normalEstimation(pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud, const float &search_radius){
        if(point_cloud->points.size() == 0){
            return nullptr;
        }
        else{
            pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(new pcl::PointCloud<pcl::Normal>);
            pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;
            normalEstimation.setInputCloud(point_cloud);

            pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
            normalEstimation.setSearchMethod(tree);

            pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud< pcl::Normal>);
            normalEstimation.setRadiusSearch (search_radius);
            normalEstimation.compute (*normals);

            pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(new pcl::PointCloud<pcl::PointNormal>);
            pcl::concatenateFields(*point_cloud, *normals, *cloud_with_normals);

            return cloud_with_normals;
        }
    }

    void pointCloud::computeShotFeatures(pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals)
    {
       if(cloud_with_normals->points.size() == 0)
       {
        std::exit(EXIT_FAILURE);
       } 
       else{
        pcl::SHOTEstimationOMP<pcl::PointNormal, pcl::PointNormal, pcl::SHOT352> shot;
        shot.setInputCloud(cloud_with_normals);
        shot.setInputNormals(cloud_with_normals);
        pcl::search::KdTree<pcl::PointNormal>::Ptr tree(new pcl::search::KdTree<pcl::PointNormal>());
        shot.setSearchMethod(tree);

        //Estimate SHOT descriptors
        pcl::PointCloud<pcl::SHOT352>::Ptr shot_descriptor(new pcl::PointCloud<pcl::SHOT352>);
        shot.compute(*shot_descriptor);

        if(cloud_with_normals->size() == shot_descriptor->size())
        {
            for(uint32_t i = 0;i < cloud_with_normals->points.size();++i)
            {
                Eigen::Vector3d point, normal;
                std::array<double, 352> desc{};
                point[0] = cloud_with_normals->points.at(i).x;
                point[1] = cloud_with_normals->points.at(i).y;
                point[2] = cloud_with_normals->points.at(i).z;
                
                normal[0] = cloud_with_normals->points.at(i).normal_x;
                normal[1] = cloud_with_normals->points.at(i).normal_y;
                normal[2] = cloud_with_normals->points.at(i).normal_z;

                for(uint32_t j = 0;j < shot_descriptor->at(i).descriptorSize();++j)
                {
                    desc[j] = shot_descriptor->at(i).descriptor[j];
                }

                m_point_cloud.emplace_back(Point(point, normal, desc));
            }
        }
       }
    }
}