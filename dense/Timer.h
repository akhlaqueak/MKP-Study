/******************************************************************************
 * timer.h 
 * record the elapsed time in microseconds (10^{-6} second)
 *****************************************************************************/

#ifndef _TIMER_LJ_
#define _TIMER_LJ_

#include <cstdlib>
#include <sys/time.h>
#include <chrono>
using namespace std;
class Timer
{
#define TIME_NOW chrono::steady_clock::now()

public:
	Timer() : m_start(TIME_NOW) {}
	void restart() { m_start = TIME_NOW; }
	long long elapsed()
	{
		return chrono::duration_cast<chrono::microseconds>(TIME_NOW - m_start).count();
	}
	void tick(){
		tic = TIME_NOW;
	}
	void tock(){
		toc+=chrono::duration_cast<chrono::nanoseconds>(TIME_NOW - tic).count();
	}
	double ticktock(){
		return toc/1e9;
	}
private:
	std::chrono::steady_clock::time_point m_start, tic;
	unsigned long long toc=0;
};

#endif
