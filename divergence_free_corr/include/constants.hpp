#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include<thread>

const double SEARCH_RADIUS = 0.003;
const uint32_t NUMBER_OF_SAMPLE_POINTS = 400;
constexpr double SIGMA = 0.01;
const uint32_t MAX_NUMBER_OF_VELOCITY_BASIS = 50;
const uint32_t MAX_NUMBER_OF_ITERS = 100;
const uint32_t NUMBER_OF_TIME_STEPS = 1;
const double TOLERANCE = 1e-06;
constexpr double r0 = 0.01;
const uint32_t MAX_NUMBER_OF_THREADS = std::thread::hardware_concurrency();

#endif //CONSTANTS_HPP