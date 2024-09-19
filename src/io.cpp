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

    bool isValidDescriptor(const pcl::SHOT352& descriptor) {
        for (const auto& value : descriptor.descriptor) {
            if (std::isnan(value) || std::isinf(value)) {
                return false;
            }
        }
        return true;
    }

    std::unique_ptr<std::vector<std::pair<adi::Point, adi::Point>>> computeCorrespondences(
    std::vector<adi::Point> source_point_cloud,
    std::vector<adi::Point> target_point_cloud)
    {
        assert(source_point_cloud.size() == target_point_cloud.size());
        
        std::unique_ptr<std::vector<std::pair<adi::Point, adi::Point>>> correspondences = std::make_unique<std::vector<std::pair<adi::Point, adi::Point>>>();
        
        // Convert target point cloud to PCL cloud with SHOT descriptors
        pcl::PointCloud<pcl::SHOT352>::Ptr target_descriptors(new pcl::PointCloud<pcl::SHOT352>);
        for (adi::Point& point : target_point_cloud)
        {
            target_descriptors->push_back(point.convertToPCLDescriptor());
        }

        // Create a KdTree for descriptor matching
        pcl::KdTreeFLANN<pcl::SHOT352> kdtree;
        kdtree.setInputCloud(target_descriptors);

        // Find correspondences
        for (adi::Point& source_point : source_point_cloud)
        {
            pcl::SHOT352 source_descriptor = source_point.convertToPCLDescriptor();
            if (!isValidDescriptor(source_descriptor)) {
            std::cerr << "Invalid source descriptor." << std::endl;
            continue;
            }
            std::vector<int> nn_indices(1);
            std::vector<float> nn_distances(1);
            
            // Perform nearest neighbor search
            if (kdtree.nearestKSearch(source_descriptor, 1, nn_indices, nn_distances) > 0)
            {
                int target_index = nn_indices[0];
                correspondences->push_back(std::make_pair(source_point, target_point_cloud[target_index]));
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
        std::cout << "Size of pointcloud after loading: " << point_cloud->size() << std::endl;
        this->computeShotFeatures(point_cloud, search_radius);
        std::cout << "SHOT features computed" << std::endl;
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
        // Perform normal estimation
        pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals = this->normalEstimation(cloud, search_radius);
        
        std::cout << "Size of cloud with normals: " << cloud_with_normals->size() << std::endl; 

        if (cloud_with_normals->empty())
        {
            std::cerr << "Error: Cloud with normals is empty. Exiting." << std::endl;
            std::exit(EXIT_FAILURE);
        } 
        else
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

            // Check if sizes match
            if (cloud_with_normals->size() != shot_descriptor->size())
            {
                std::cerr << "Error: Sizes of cloud with normals and SHOT descriptors do not match." << std::endl;
                std::exit(EXIT_FAILURE);
            }

            // Store features
            for (std::size_t i = 0; i < cloud_with_normals->points.size(); ++i)
            {
                Eigen::Vector3d point, normal;
                std::array<double, 352> desc{};
                point[0] = cloud_with_normals->points[i].x;
                point[1] = cloud_with_normals->points[i].y;
                point[2] = cloud_with_normals->points[i].z;

                normal[0] = cloud_with_normals->points[i].normal_x;
                normal[1] = cloud_with_normals->points[i].normal_y;
                normal[2] = cloud_with_normals->points[i].normal_z;

                // Copy descriptor values
                for (std::size_t j = 0; j < shot_descriptor->at(i).descriptorSize(); ++j)
                {
                    desc[j] = shot_descriptor->at(i).descriptor[j];
                }

                // Store the point, normal, and descriptor in your data structure
                m_point_cloud.emplace_back(point, normal, desc);
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

    std::unique_ptr<std::vector<adi::Point>> pointCloud::samplePointCloud(const uint32_t max_number_of_points)
    {   
        if (max_number_of_points == 0) {
            return std::make_unique<std::vector<adi::Point>>();
        }

        // Create a unique pointer for the subsampled point cloud
        std::vector<adi::Point> subsampled_point_cloud;
        subsampled_point_cloud.reserve(max_number_of_points);

        if (m_point_cloud.size() >= max_number_of_points) {
            // Initialize variables
            std::vector<bool> is_sampled(m_point_cloud.size(), false);
            std::vector<double> min_distances(m_point_cloud.size(), std::numeric_limits<double>::max());

            // Randomly pick the first point
            int random_index = std::rand() % m_point_cloud.size();
            subsampled_point_cloud.push_back(m_point_cloud[random_index]);
            is_sampled[random_index] = true;

            // Update distances for the first point
            std::priority_queue<std::pair<double, size_t>, std::vector<std::pair<double, size_t>>, FarthestPointComparator> max_heap;
            for (size_t j = 0; j < m_point_cloud.size(); ++j) {
            if (j != random_index) {
                min_distances[j] = (m_point_cloud[j].s_point - m_point_cloud[random_index].s_point).norm();
                max_heap.emplace(min_distances[j], j);
            }
            }

            // Iteratively sample the farthest points
            for (size_t i = 1; i < max_number_of_points; ++i) {
            // Find the farthest point using the max-heap
            size_t farthest_index;
            do {
                farthest_index = max_heap.top().second;
                max_heap.pop();
            } while (is_sampled[farthest_index]);  // Ensure we pick an unsampled point

            subsampled_point_cloud.push_back(m_point_cloud[farthest_index]);
            is_sampled[farthest_index] = true;

            // Update the minimum distances and max-heap for remaining points
            const Eigen::Vector3d &new_sampled_point = m_point_cloud[farthest_index].s_point;
            for (size_t j = 0; j < m_point_cloud.size(); ++j) {
                if (!is_sampled[j]) {
                double distance = (m_point_cloud[j].s_point - new_sampled_point).norm();
                if (distance < min_distances[j]) {
                    min_distances[j] = distance;
                    max_heap.emplace(min_distances[j], j);  // Push updated distance into heap
                }
                }
            }
            }

            return std::make_unique<std::vector<Point>>(std::move(subsampled_point_cloud));
        } else {
            return std::make_unique<std::vector<Point>>(std::move(m_point_cloud));
        }
    }
}
