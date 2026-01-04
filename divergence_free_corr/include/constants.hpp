#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <thread>

const double SEARCH_RADIUS =
    0.0576; // This value is to be adjusted accoring to the sampled points size
constexpr int NUMBER_OF_SAMPLE_POINTS = 500;
constexpr double SIGMA = 0.1;
constexpr int MAX_NUMBER_OF_VELOCITY_BASIS = 50;
const uint32_t MAX_NUMBER_OF_ITERS = 1; // 100
const uint32_t NUMBER_OF_TIME_STEPS = 1000;
const double TOLERANCE = 1e-06;
constexpr double r0 = 0.01;
const uint32_t MAX_NUMBER_OF_THREADS = std::thread::hardware_concurrency();
constexpr int TENSOR_DEPTH = 9;

#endif // CONSTANTS_HPP
