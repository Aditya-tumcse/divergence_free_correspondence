#include "io.hpp"

#include<cassert>

namespace utilities{
    double eucledianDistance(const Eigen::Vector3d &point_1, const Eigen::Vector3d &point_2)
    {
        return (std::sqrt(std::pow((point_1[0] - point_2[0]), 2) + std::pow((point_1[1] - point_2[1]), 2) + std::pow((point_1[2] - point_2[2]), 2)));
    }

    Eigen::MatrixXd computeDistanceMatrix(const std::vector<Eigen::Vector3d> points){
        Eigen::MatrixXd distance_matrix = Eigen::MatrixXd::Zero(points.size(), points.size());
        for(int i = 0;i < points.size();++i){
            for(int j = 0;j < points.size();++j){
                double eucledian_distance = eucledianDistance(points[i], points[j]);
                distance_matrix(i,j) = eucledian_distance;
            }
        }
        return distance_matrix;
    }

    double computeL2Norm(const std::array<double, 352> &descriptor_1, const std::array<double, 352>  &descriptor_2)
    {   
        std::array<double, 352> difference{};
        for(uint32_t i = 0;i < descriptor_1.size();++i)
        {   
            difference[i] = descriptor_1[i] - descriptor_2[i];
        }
        double sum = 0.0;
        for(double element : difference)
        {
            sum += element * element;
        }
        return std::sqrt(sum);
    }

    std::vector<std::pair<adi::Point, adi::Point>> computeCorrespondences(std::vector<adi::Point> source_point_cloud,std::vector<adi::Point> target_point_cloud)
    {

        assert(source_point_cloud.size() == target_point_cloud.size());
        std::vector<std::pair<adi::Point, adi::Point>> correspondences;
        pcl::PointCloud<pcl::SHOT352>::Ptr target_descriptors(new pcl::PointCloud<pcl::SHOT352>);
        for(adi::Point &point : target_point_cloud)
        {
            target_descriptors->push_back(point.convertToPCLDescriptor());
        }
        pcl::KdTreeFLANN<pcl::SHOT352> kdtree;
        kdtree.setInputCloud(target_descriptors);

        // Find correspondence
        for(auto &source_point : source_point_cloud)
        {
            pcl::SHOT352 source_descriptor = source_point.convertToPCLDescriptor();
            std::vector<int> nn_indices(1);
            std::vector<float> nn_distances(1);

            if(kdtree.nearestKSearch(source_descriptor, 1, nn_indices, nn_distances) > 0)
            {
                int target_index = nn_indices[0];
                correspondences.emplace_back(std::make_pair(source_point, target_point_cloud[target_index]));
            }
        }
        return correspondences;
    }
}


namespace adi{

    pcl::SHOT352 Point::convertToPCLDescriptor()
     {
        pcl::SHOT352 descriptor;
        for(uint32_t i=0;i < s_descriptor.size();++i)
        {
            descriptor.descriptor[i] = static_cast<float>(s_descriptor[i]);
        }
        return descriptor;
    }

    pointCloud::pointCloud(const std::string &point_cloud_path, const double &search_radius){
        pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud = this->readPointCloud(point_cloud_path);
        this->computeShotFeatures(point_cloud, search_radius);
    }

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

    pcl::PointCloud<pcl::PointNormal>::Ptr pointCloud::normalEstimation(pcl::PointCloud<pcl::PointXYZ>::Ptr point_cloud, const double &search_radius){
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

    void pointCloud::computeShotFeatures(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, const double &search_radius)
    {
       pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals = this->normalEstimation(cloud, search_radius); 
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

    const std::vector<Eigen::Vector3d> pointCloud::extractPoints(std::vector<adi::Point> point_cl)
    {
        std::vector<Eigen::Vector3d> points;
        points.reserve(point_cl.size());

        for(int i = 0;i < point_cl.size();++i)
        {
            points.emplace_back(point_cl[i].s_point);
        }
        return points;
    }

    uint32_t pointCloud::getColumnIdOfTheFarthestSample(const Eigen::RowVectorXd &row)
    {
        std::vector<double> vec(row.data(), row.data() + row.size());
        auto maxIndex = [](const std::vector<double> &vec){
            auto it = std::max_element(vec.begin(), vec.end());
            return static_cast<uint32_t>(std::distance(vec.begin(), it));
        };

        return maxIndex(vec);
    }

    std::vector<adi::Point> pointCloud::samplePointCloud(const uint32_t max_number_of_points)
    {   
        assert(max_number_of_points > 1);
        Eigen::MatrixXd distance_matrix = utilities::computeDistanceMatrix(this->extractPoints(m_point_cloud));
        std::vector<adi::Point> subsampled_point_cloud;
        subsampled_point_cloud.reserve(max_number_of_points);

        std::vector<adi::Point> og_point_cloud_copy = m_point_cloud;
        if(m_point_cloud.size() != 0){
            int random_index = std::rand() % m_point_cloud.size();
            adi::Point anchor_point = og_point_cloud_copy[random_index];
            subsampled_point_cloud.emplace_back(anchor_point);
            for(int i = 0;i < max_number_of_points - 1;++i)
            {
                const int row_id = random_index;
                const double farthest_point_pt = distance_matrix.row(row_id).maxCoeff();
                const int index_of_farthest_point = this->getColumnIdOfTheFarthestSample(distance_matrix.row(row_id));
                const adi::Point farthest_point = og_point_cloud_copy[index_of_farthest_point];
                subsampled_point_cloud.emplace_back(farthest_point);
                distance_matrix.row(row_id).setConstant(-100);
                distance_matrix.col(row_id).setConstant(-100);
                std::swap(og_point_cloud_copy[index_of_farthest_point], og_point_cloud_copy.back());
                og_point_cloud_copy.pop_back();
                random_index = index_of_farthest_point;
            }
        }
        else{
            std::cout << "Input point cloud is empty" << std::endl;
        }
        return subsampled_point_cloud;
    }
}