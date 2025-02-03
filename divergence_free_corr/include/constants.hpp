#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <thread>

constexpr double SEARCH_RADIUS =
    0.0576; // This value is to be adjusted accoring to the sampled points size
constexpr int NUMBER_OF_SAMPLE_POINTS = 500;
constexpr double SIGMA = 0.1;
constexpr int MAX_NUMBER_OF_VELOCITY_BASIS = 50;
constexpr uint32_t MAX_NUMBER_OF_ITERS = 1; // 100
constexpr uint32_t NUMBER_OF_TIME_STEPS = 1000;
constexpr double TOLERANCE = 1e-06;
constexpr int TENSOR_DEPTH = 9;
constexpr double THRESHOLD_CORR = 0.01;
const size_t MAX_NUMBER_OF_THREADS = std::thread::hardware_concurrency();

#endif // CONSTANTS_HPP