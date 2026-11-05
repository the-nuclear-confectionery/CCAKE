#ifndef UTILITIES_H
#define UTILITIES_H

#include <cstdlib>
#include <cstdio>
#include <string>
#ifdef __linux__
#include <sys/sysinfo.h>
#include <unistd.h>       // sysconf, _SC_PAGESIZE
#endif

/// Lightweight heap probe: reads VmRSS from /proc/self/statm (one syscall, ~2 µs).
/// Returns RSS in MB.  Use the HEAP_PROBE macro for annotated logging.
inline double get_rss_mb() {
#ifdef __linux__
    long pages = 0;
    FILE* f = std::fopen("/proc/self/statm", "r");
    if (f) { long dummy; std::fscanf(f, "%ld %ld", &dummy, &pages); std::fclose(f); }
    return pages * (sysconf(_SC_PAGESIZE) / 1048576.0);
#else
    return 0.0;
#endif
}

/// Print an annotated RSS delta.  Usage:
///   double _rss0 = get_rss_mb();
///   ... code ...
///   HEAP_PROBE("label", _rss0);
#define HEAP_PROBE(label, rss_before) do { \
    double _now = get_rss_mb(); \
    std::fprintf(stderr, "[HEAP] %-42s  RSS: %8.1f MB   Δ: %+8.1f MB\n", \
                 (label), _now, _now - (rss_before)); \
} while(0)

/// Global dfinitions used in the Cabana and Kokkos libraries
#include <Cabana_Core.hpp>

#define VECTOR_LENGTH 8 ///< The length of the vector used in the Cabana library
                        ///< AMD EPYC 7763 has AVX2 (256-bit = 4 doubles); 8 = 2x unroll,
                        ///< good balance of vectorization and cache reuse for double-precision SPH.
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

// Cabana >= 0.8.0 dropped the Kokkos::Device template parameter for AoSoA/Slice in
// favour of a Kokkos memory space. Kokkos::View also accepts a bare memory space, so
// using the memory space here works in both the AoSoA and View contexts (CPU and GPU).
using DeviceType = MemorySpace;
using HostType = Kokkos::HostSpace;

/// Alias for the particle list
using ListAlgorithm = Cabana::FullNeighborTag;
using ListLayout = Cabana::VerletLayoutCSR; ///< CSR packs neighbors contiguously per particle;
                                            ///< VerletLayout2D pads every row to max_neighbors,
                                            ///< blowing up to ~N*max_neigh*4B when neighbor counts vary.
using ListType = Cabana::VerletList<MemorySpace, ListAlgorithm, ListLayout, Cabana::TeamOpTag>;

/// @brief Sorts rows 0..number_of_rows-1 of a CSR neighbour list into ascending index order.
/// @details With per-particle radii the (patched) Cabana fill adds cross-neighbours into
/// other particles' rows from racing threads, so a row's ORDER depends on thread
/// scheduling while its membership does not. Every SPH or WVT sum over a row is a
/// floating-point sum, and FP addition is not associative, so an unsorted row makes the
/// result change from run to run. Sorting makes the order canonical. Rows are short
/// (~50-200 entries), so an in-place insertion sort is used.
inline void sort_neighbor_rows(ListType& list, const int number_of_rows)
{
  auto neighbor_counts  = list._data.counts;
  auto neighbor_offsets = list._data.offsets;
  auto neighbor_ids     = list._data.neighbors;
  Kokkos::parallel_for("sort_neighbor_rows", Kokkos::RangePolicy<ExecutionSpace>(0, number_of_rows),
    KOKKOS_LAMBDA(const int i){
      const int row_begin = neighbor_offsets(i);
      const int row_length = neighbor_counts(i);
      for (int a = 1; a < row_length; ++a) {
        const int value = neighbor_ids(row_begin + a);
        int b = a - 1;
        while (b >= 0 && neighbor_ids(row_begin + b) > value) {
          neighbor_ids(row_begin + b + 1) = neighbor_ids(row_begin + b);
          --b;
        }
        neighbor_ids(row_begin + b + 1) = value;
      }
    });
  Kokkos::fence();
}

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
    (*Ainv)[0][0] = 1/(*A)[0][0];
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