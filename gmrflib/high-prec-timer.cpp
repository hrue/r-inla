//#include <bits/stdc++.h>
#include <chrono>

using namespace std;
extern "C" double GMRFLib_timer_chrono(void);

double GMRFLib_timer_chrono(void)
{
	static auto start = chrono::high_resolution_clock::now();
	auto now = chrono::high_resolution_clock::now();
	return ((chrono::duration_cast<chrono::nanoseconds>(now - start).count()) * 1.0E-09);
}
