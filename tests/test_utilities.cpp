#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "utilities.hpp"

using Catch::Matchers::WithinRel;

TEST_CASE("eucledianDistance: identical points return zero", "[utilities]") {
    Eigen::Vector3d p{1.0, 2.0, 3.0};
    REQUIRE_THAT(utilities::eucledianDistance(p, p), WithinRel(0.0, 1e-12));
}

TEST_CASE("eucledianDistance: unit axis distances", "[utilities]") {
    Eigen::Vector3d origin{0.0, 0.0, 0.0};
    REQUIRE_THAT(utilities::eucledianDistance(origin, {1.0, 0.0, 0.0}), WithinRel(1.0, 1e-12));
    REQUIRE_THAT(utilities::eucledianDistance(origin, {0.0, 1.0, 0.0}), WithinRel(1.0, 1e-12));
    REQUIRE_THAT(utilities::eucledianDistance(origin, {0.0, 0.0, 1.0}), WithinRel(1.0, 1e-12));
}

TEST_CASE("eucledianDistance: 3-4-5 right triangle", "[utilities]") {
    Eigen::Vector3d a{0.0, 0.0, 0.0};
    Eigen::Vector3d b{3.0, 4.0, 0.0};
    REQUIRE_THAT(utilities::eucledianDistance(a, b), WithinRel(5.0, 1e-12));
}

TEST_CASE("eucledianDistance: symmetric", "[utilities]") {
    Eigen::Vector3d a{1.0, 2.0, 3.0};
    Eigen::Vector3d b{4.0, 6.0, 8.0};
    REQUIRE_THAT(utilities::eucledianDistance(a, b),
                 WithinRel(utilities::eucledianDistance(b, a), 1e-12));
}

TEST_CASE("computeL2Norm: identical descriptors return zero", "[utilities]") {
    std::array<double, 352> d{};
    d.fill(1.0);
    REQUIRE_THAT(utilities::computeL2Norm(d, d), WithinRel(0.0, 1e-12));
}

TEST_CASE("computeL2Norm: single nonzero element", "[utilities]") {
    std::array<double, 352> a{};
    std::array<double, 352> b{};
    a[0] = 3.0;
    b[0] = 0.0;
    REQUIRE_THAT(utilities::computeL2Norm(a, b), WithinRel(3.0, 1e-12));
}

TEST_CASE("computeL2Norm: symmetric", "[utilities]") {
    std::array<double, 352> a{};
    std::array<double, 352> b{};
    a[10] = 1.0;
    b[10] = 4.0;
    a[20] = 2.0;
    b[20] = 6.0;
    REQUIRE_THAT(utilities::computeL2Norm(a, b),
                 WithinRel(utilities::computeL2Norm(b, a), 1e-12));
}

TEST_CASE("toEigenMatrix / toPointCloud roundtrip", "[utilities]") {
    std::vector<adi::Point> original;
    original.emplace_back(Eigen::Vector3d{1.0, 2.0, 3.0}, Eigen::Vector3d{0, 0, 1},
                          std::array<double, 352>{});
    original.emplace_back(Eigen::Vector3d{4.0, 5.0, 6.0}, Eigen::Vector3d{0, 1, 0},
                          std::array<double, 352>{});

    Eigen::MatrixXd matrix = utilities::toEigenMatrix(original);
    REQUIRE(matrix.rows() == 2);
    REQUIRE(matrix.cols() == 3);

    std::vector<adi::Point> recovered = utilities::toPointCloud(matrix);
    REQUIRE(recovered.size() == 2);
    REQUIRE_THAT(recovered[0].s_point.x(), WithinRel(1.0, 1e-12));
    REQUIRE_THAT(recovered[1].s_point.z(), WithinRel(6.0, 1e-12));
}

TEST_CASE("toPointCloud: NaN rows are silently dropped", "[utilities]") {
    Eigen::MatrixXd m(3, 3);
    m << 1.0, 2.0, 3.0,
         std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0,
         4.0, 5.0, 6.0;

    auto cloud = utilities::toPointCloud(m);
    REQUIRE(cloud.size() == 2);
}
