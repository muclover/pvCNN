#pragma once

#include <assert.h>
// #include <iostream>
#include <chrono>

class Timer
{
public:
    Timer() : total_time(0), status(false)
    {
        // start();
    }
    void start();
    void stop();
    void clear()
    {
        total_time = 0;
        status = false;
    }
    double elapse_time()
    {
        assert(status == false);
        return total_time;
    }

    // ~Timer()
    // {
    // std::cout << "time : " << total_time << std::endl;
    // stop();
    // }

private:
    std::chrono::high_resolution_clock::time_point t;
    double total_time;
    bool status;
};