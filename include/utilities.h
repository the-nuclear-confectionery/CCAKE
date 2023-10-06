#ifndef UTILITIES_H
#define UTILITIES_H

#include <cstdlib>
#ifdef __linux__
#include <sys/sysinfo.h>
#endif

class Utilities
{
public:
    static unsigned long get_free_memory();
    static unsigned long get_total_memory();
    static void lack_of_memory_stop();
};
#endif // UTILITIES_H