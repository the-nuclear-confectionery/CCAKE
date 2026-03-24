#ifndef SYSTEM_STATE_INT_H
#define SYSTEM_STATE_INT_H

#include <string>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Cabana_LinkedCellList.hpp>

#include "kernel.h"

class SystemState {
public:
    struct Particle {
        int    id;
        double t;

        double r[3];

        double p, T, muB, muS, muQ;
        double e, rhoB, rhoS, rhoQ, s;

        double eta_pi, zeta_Pi, tau_Pi, tau_pi;
        double theta;
        double invRe_shear, invRe_bulk;
        double kn_shear, kn_bulk;

        double shv[9];   // 00 11 22 01 02 12 13 23 33

        double u[3];
        double gamma;

        int Freeze;
        std::string eos;
        double causality;

        double sph_mass;
    };

    // Compact struct holding only the fields needed for SPH interpolation.
    // Stored in cell-sorted order to eliminate permutation indirection and
    // improve cache utilisation (80 B vs ~280 B for the full Particle).
    struct SphFields {
        double r[3];
        double e, s, rhoB, rhoS, rhoQ, T;
        double sph_mass;
    };

    struct InterpResult {
        double sigma_star = 0.0;
        double e     = 0.0;
        double s     = 0.0;
        double rhoB  = 0.0;
        double rhoS  = 0.0;
        double rhoQ  = 0.0;
        double T     = 0.0;
    };

    // Convenience type for batch queries
    struct Query {
        double t;
        double x[3];
        double hT;
    };

    using CellList = Cabana::LinkedCellList<Kokkos::HostSpace>;

    struct TimeStep {
        double time;
        std::vector<Particle>  particles;  // full data (kept for non-SPH access)
        std::vector<SphFields> sph;        // cell-sorted compact SPH data

        // built once at load time for fast spatial queries
        CellList cell_list;
        double   lo[3];       // grid origin (with padding)
        double   cell_size;   // cell size used to build the list
    };

public:
    // configuration
    int    D;
    double t0;
    double dt;
    int    output_steps;
    std::string output_dir;

    // single-snapshot data (legacy load() interface)
    double file_time;
    std::vector<Particle> particles;

    // full evolution (populated by load_all())
    std::vector<TimeStep> evolution;

public:
    SystemState(int D_from_config,
                const std::string& output_dir,
                double t0, double dt,
                int output_steps);

    // load single snapshot at time t (legacy)
    void load(double t);

    // load all snapshots from output_dir; h_cell is the cell size for the
    // linked cell list (should be >= 2*hT used at query time for best coverage)
    void load_all(double h_cell = 0.5);

    // SPH interpolation using this->particles (legacy, no cell list)
    InterpResult interpolate(const double* x, double hT) const;

    // SPH in space + linear in time interpolation (sequential, uses bracket cache)
    InterpResult interpolate_at(double t, const double* x, double hT) const;

    // Parallel batch version — thread-safe, no shared state
    std::vector<InterpResult> interpolate_many(
        const std::vector<Query>& queries) const;

private:
    // Bracket cache for sequential interpolate_at calls (not thread-safe)
    mutable std::size_t last_lo_idx_ = 0;

    std::string filename_for_time(double t) const;
    std::string filename_for_index(long long idx) const;
    void load_into(const std::string& filename, TimeStep& ts,
                   double h_cell) const;

    // Core implementation: find bracket by index and interpolate (no cache)
    InterpResult interpolate_at_impl(double t, const double* x, double hT,
                                     std::size_t lo_idx) const;

    template<unsigned int DD>
    InterpResult interpolate_impl(const TimeStep& ts,
                                  const double* x, double hT) const;

    // fallback for single-snapshot (legacy): O(N) loop, no cell list
    template<unsigned int DD>
    InterpResult interpolate_impl_slow(const std::vector<Particle>& parts,
                                       const double* x, double hT) const;
};

#endif
