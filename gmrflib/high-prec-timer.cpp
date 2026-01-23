//#include <bits/stdc++.h>
#include <chrono>

using namespace std;
extern "C" {
	double GMRFLib_timer_chrono(void);
	double GMRFLib_timer_chrono_mono(void);
}

double GMRFLib_timer_chrono(void)
{
	static auto start = chrono::high_resolution_clock::now();
	auto now = chrono::high_resolution_clock::now();
	return ((chrono::duration_cast<chrono::nanoseconds>(now - start).count()) * 1.0E-09);
}

double GMRFLib_timer_chrono_mono(void)
{
    static auto start = std::chrono::steady_clock::now();
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(now - start).count() * 1.0E-9;
}
