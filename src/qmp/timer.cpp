#include "timer.hpp"

void Timer::start()
{
    assert(status == false);
    t = std::chrono::high_resolution_clock::now();
    status = true;
}

void Timer::stop()
{
    assert(status == true);
    auto t_end = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t);
    total_time += time.count();
    status = false;
}