#ifndef UTILITIES_H
#define UTILITIES_H

#include <cstdlib>
#ifdef __linux__
#include <sys/sysinfo.h>
#endif

/// Global dfinitions used in the Cabana and Kokkos libraries
#include <Cabana_Core.hpp>

#define VECTOR_LENGTH 32 ///< The length of the vector used in the Cabana library
/// Choose between CPU or GPU execution
#ifdef __CUDACC__
using MemorySpace = Kokkos::CudaSpace;
using ExecutionSpace = Kokkos::Cuda;
#else
using MemorySpace = Kokkos::HostSpace;
using ExecutionSpace = Kokkos::OpenMP;
#endif

using DeviceType = Kokkos::Device<ExecutionSpace, MemorySpace>;

/// Alias for the particle list
using ListAlgorithm = Cabana::FullNeighborTag;
using ListLayout = Cabana::VerletLayoutCSR; ///< There are other options for the list layout. //TODO: Check which one is the best
using ListType = Cabana::VerletList<MemorySpace, ListAlgorithm, ListLayout>;

class Utilities
{
public:
    static unsigned long get_free_memory();
    static unsigned long get_total_memory();
    static void lack_of_memory_stop();
};
#endif // UTILITIES_H