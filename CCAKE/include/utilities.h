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
    #ifdef DEBUG
        using MemorySpace = Kokkos::CudaUVMSpace; //Use managed space if debugging the code - Should we keep this for final release?
    #else
        using MemorySpace = Kokkos::CudaUVMSpace;
    #endif
using ExecutionSpace = Kokkos::Cuda;
#else
using MemorySpace = Kokkos::HostSpace;
using ExecutionSpace = Kokkos::OpenMP;
#endif

using DeviceType = Kokkos::Device<ExecutionSpace, MemorySpace>;
using HostType = Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>;

/// Alias for the particle list
using ListAlgorithm = Cabana::FullNeighborTag;
using ListLayout = Cabana::VerletLayout2D; ///< There are other options for the list layout. //TODO: Check which one is the best
using ListType = Cabana::VerletList<MemorySpace, ListAlgorithm, ListLayout,Cabana::TeamOpTag>;

namespace Utilities
{
    static unsigned long get_free_memory();
    static unsigned long get_total_memory();
    static void lack_of_memory_stop();

    template<unsigned int D>
    KOKKOS_INLINE_FUNCTION
    void inverse(double (*A)[D][D], double (*Ainv)[D][D]);
};

template<>
KOKKOS_INLINE_FUNCTION
void Utilities::inverse<1>(double (*A)[1][1], double (*Ainv)[1][1])
{
    (*Ainv)[0][0] = 1/(*A)[1][1];
}

template<>
KOKKOS_INLINE_FUNCTION
void Utilities::inverse<2>(double (*A)[2][2], double (*Ainv)[2][2])
{
    double det = (*A)[0][0]*(*A)[1][1] - (*A)[0][1]*(*A)[1][0];
    (*Ainv)[0][0] = (*A)[1][1]/det;
    (*Ainv)[0][1] = -(*A)[0][1]/det;
    (*Ainv)[1][0] = -(*A)[1][0]/det;
    (*Ainv)[1][1] = (*A)[0][0]/det;
}

template<>
KOKKOS_INLINE_FUNCTION
void Utilities::inverse<3>(double (*A)[3][3], double (*Ainv)[3][3])
{
    double det = (*A)[0][0]*(*A)[1][1]*(*A)[2][2] + (*A)[0][1]*(*A)[1][2]*(*A)[2][0] +
                 (*A)[0][2]*(*A)[1][0]*(*A)[2][1] - (*A)[0][2]*(*A)[1][1]*(*A)[2][0] -
                 (*A)[0][1]*(*A)[1][0]*(*A)[2][2] - (*A)[0][0]*(*A)[1][2]*(*A)[2][1];
    (*Ainv)[0][0] = ((*A)[1][1]*(*A)[2][2] - (*A)[1][2]*(*A)[2][1])/det;
    (*Ainv)[0][1] = ((*A)[0][2]*(*A)[2][1] - (*A)[0][1]*(*A)[2][2])/det;
    (*Ainv)[0][2] = ((*A)[0][1]*(*A)[1][2] - (*A)[0][2]*(*A)[1][1])/det;
    (*Ainv)[1][0] = ((*A)[1][2]*(*A)[2][0] - (*A)[1][0]*(*A)[2][2])/det;
    (*Ainv)[1][1] = ((*A)[0][0]*(*A)[2][2] - (*A)[0][2]*(*A)[2][0])/det;
    (*Ainv)[1][2] = ((*A)[0][2]*(*A)[1][0] - (*A)[0][0]*(*A)[1][2])/det;
    (*Ainv)[2][0] = ((*A)[1][0]*(*A)[2][1] - (*A)[1][1]*(*A)[2][0])/det;
    (*Ainv)[2][1] = ((*A)[0][1]*(*A)[2][0] - (*A)[0][0]*(*A)[2][1])/det;
    (*Ainv)[2][2] = ((*A)[0][0]*(*A)[1][1] - (*A)[0][1]*(*A)[1][0])/det;
}
#endif // UTILITIES_H