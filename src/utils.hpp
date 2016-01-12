/*!
 * @file utils.hpp
 *
 * @brief utils header file
 */

#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <memory>
#include <future>
#include <sys/time.h>

class Timer {
public:

    /*!
     * @brief Timer constructor
     */
    Timer();

    /*!
     * @brief Method for starting the timer
     * @details Method records the current time into timeval_ and unpauses
     * the timer.
     */
    void start();

    /*!
     * @brief Method for stopping the timer
     * @details Method subtracts the current time from the time in timeval_
     * and adds the difference to time_ if the timer is not paused.
     */
    void stop();

    /*!
     * @brief Method for reseting the timer
     * @details Method resets the variable time_ and starts the timer.
     */
    void reset();

    /*!
     * @brief Method for elapse time printing
     * @details Method prints to stderr the elapsed time in seconds in following
     * format: [location][message] time (s).
     *
     * @param [in] location used for location (function/module) the timer is in
     * @param [in] message used for action the timer was timing
     */
    void print(const char* location, const char* message) const;

private:

    bool paused_;
    uint64_t time_;
    timeval timeval_;
};
