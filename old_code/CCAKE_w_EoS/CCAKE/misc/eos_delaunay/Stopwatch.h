#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <ctime>

class Stopwatch
{
	private:
		time_t start, end;
		time_t current_elapsed_time;
	public:
		Stopwatch() {start=clock(); end=0; current_elapsed_time=0;}
		void Start() {start=clock();}
		void Stop() {end=clock(); current_elapsed_time+=(end - start);}
		void Reset() {start=clock(); end=0; current_elapsed_time=0;}
		double printTime() {return ((double)current_elapsed_time) / CLOCKS_PER_SEC;}
};

#endif
