/*!
 * @file utils.cpp
 *
 * @brief utils source file
 */

#include <assert.h>

#include "utils.hpp"

Timer::Timer()
    : paused_(false), time_(0), timeval_() {
}

void Timer::start() {

    gettimeofday(&timeval_, nullptr);
    paused_ = false;
}

void Timer::stop() {

    if (paused_) return;

    timeval stop;
    gettimeofday(&stop, nullptr);

    time_ += ((stop.tv_sec - timeval_.tv_sec) * 1000000L + stop.tv_usec)
        - timeval_.tv_usec;
    paused_ = true;
}

void Timer::reset() {

    gettimeofday(&timeval_, nullptr);
    time_ = 0;
    paused_ = false;
}

void Timer::print(const char* location, const char* message) const {

    fprintf(stderr, "[%s][%s]: %.5lf s\n", location, message,
        time_ / (double) 1000000);
}
