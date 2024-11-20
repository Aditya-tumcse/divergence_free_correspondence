#ifndef RUN_HPP
#define RUN_HPP

#include "constants.hpp"
#include "deformation_field.hpp"
#include "io.hpp"
#include "matching.hpp"
#include "numerics.hpp"
#include "utilities.hpp"

/**
 * @brief Main function that runs the deformation model
 *
 * @param source_cloud
 * @param target_cloud
 */
void run(adi::pointCloud *source_cloud, adi::pointCloud *target_cloud);

#endif