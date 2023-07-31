#include "utilities.h"

///\brief Queries the amount of free memory in the system
///\return The amount of free memory (RAM) in the computer in bytes as an
///unsigned long
unsigned long Utilities::get_free_memory()
{
#ifdef __linux__
    struct sysinfo info;
    sysinfo(&info);
    return info.freeram;
#else
    std::cerr << "get_free_memory() not implemented for this OS" << std::endl;
    return 999999999999999;
#endif
}

///\brief Queries the total amount of memory in the syste,
///\return The amount of total memory intalled in the computer in bytes as an
///unsigned long
unsigned long Utilities::get_total_memory()
{
    #ifdef __linux__
    struct sysinfo info;
    sysinfo(&info);
    return info.totalram - info.freeram;
#else
    std::cerr << "get_used_memory() not implemented for this OS" << std::endl;
    return 999999999999999;
#endif
}

///\brief Stops the program if the amount of free memory is less than 1% of the
/// total memory
void Utilities::lack_of_memory_stop()
{
    unsigned long free_memory = get_free_memory();
    unsigned long total_memory = get_total_memory();
    if( (double) free_memory/ (double)total_memory > 0.99 ) std::abort();
}
