#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <chrono>
#include <ctime>
#include <ratio>

using namespace std::chrono;

class Stopwatch
{
	private:
		high_resolution_clock::time_point start, end;
		duration<double> current_elapsed_time;
	public:
		//Stopwatch() {start=clock(); end=0; current_elapsed_time=0;}
		//void Start() {start=clock();}
		//void Stop() {end=clock(); current_elapsed_time+=(end - start);}
		//void Reset() {start=clock(); end=0; current_elapsed_time=0;}
		//double printTime() {return ((double)current_elapsed_time) / CLOCKS_PER_SEC;}
		
		Stopwatch()        { start                 = high_resolution_clock::now();
			             end                   = start;
			             current_elapsed_time  = duration<double>(0);                          }

		void Start()       { start                 = high_resolution_clock::now();                 }
		void Stop()        { end                   = high_resolution_clock::now();
			             current_elapsed_time += duration_cast<duration<double>>(end - start); }
		void Reset()       { start                 = high_resolution_clock::now();
                                     end                   = start;
                                     current_elapsed_time  = duration<double>(0);                          } 
		double printTime() { return current_elapsed_time.count();                                  }
};

#endif
