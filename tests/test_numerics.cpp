#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "numerics.hpp"

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("computeResidual: zero distance returns zero", "[numerics]") {
    double src[3] = {1.0, 2.0, 3.0};
    double tgt[3] = {1.0, 2.0, 3.0};
    double result = adi::numerics::computeResidual<double>(src, tgt, 1.0);
    REQUIRE_THAT(result, WithinAbs(0.0, 1e-12));
}

TEST_CASE("computeResidual: unit distance with unit weight", "[numerics]") {
    double src[3] = {0.0, 0.0, 0.0};
    double tgt[3] = {1.0, 0.0, 0.0};
    double result = adi::numerics::computeResidual<double>(src, tgt, 1.0);
    REQUIRE_THAT(result, WithinRel(1.0, 1e-12));
}

TEST_CASE("computeResidual: scales linearly with weight", "[numerics]") {
    double src[3] = {0.0, 0.0, 0.0};
    double tgt[3] = {3.0, 4.0, 0.0};
    double r1 = adi::numerics::computeResidual<double>(src, tgt, 1.0);
    double r2 = adi::numerics::computeResidual<double>(src, tgt, 2.0);
    REQUIRE_THAT(r1, WithinRel(5.0, 1e-12));
    REQUIRE_THAT(r2, WithinRel(10.0, 1e-12));
}

TEST_CASE("GenerateBasisIndices: produces expected count", "[numerics]") {
    auto indices = adi::numerics::GenerateBasisIndices(MAX_NUMBER_OF_VELOCITY_BASIS);
    REQUIRE(indices.size() == MAX_NUMBER_OF_VELOCITY_BASIS);
}

TEST_CASE("GenerateBasisIndices: all eigenvalues are positive", "[numerics]") {
    auto indices = adi::numerics::GenerateBasisIndices(MAX_NUMBER_OF_VELOCITY_BASIS);
    for (const auto &idx : indices) {
        REQUIRE(idx.eigen_value > 0.0);
    }
}

TEST_CASE("GenerateBasisIndices: all index components are non-zero", "[numerics]") {
    auto indices = adi::numerics::GenerateBasisIndices(MAX_NUMBER_OF_VELOCITY_BASIS);
    for (const auto &idx : indices) {
        REQUIRE((idx.index_1 > 0 || idx.index_2 > 0 || idx.index_3 > 0));
    }
}
