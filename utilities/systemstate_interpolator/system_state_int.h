#ifndef SYSTEM_STATE_INT_H
#define SYSTEM_STATE_INT_H

#include <string>
#include <vector>

// use your kernel directly
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

        // ---- SPH weight (what multiplies W in sums) ----
        // If your file doesn't contain it, we default to 1.0 in load().
        double sph_mass;
    };

    // interpolation output (add/remove fields as you like)
    struct InterpResult {
        double sigma_star = 0.0;  // sum m W
        double e     = 0.0;  // sum m e W
        double s     = 0.0;  // sum m s W
        double rhoB  = 0.0;
        double rhoS  = 0.0;
        double rhoQ  = 0.0;
        double T     = 0.0;
    };

public:
    // configuration
    int    D;               // 2 or 3
    double t0;              // usually 0
    double dt;
    int    output_steps;
    std::string output_dir;

    // data
    double file_time;       // time read from file
    std::vector<Particle> particles;

public:
    SystemState(int D_from_config,
                const std::string& output_dir,
                double t0, double dt,
                int output_steps);

    // load state corresponding to time t
    void load(double t);

    // ---- SPH interpolation method (inside SystemState) ----
    // x must be length 3 (x[2] ignored if D==2). hT is smoothing length.
    InterpResult interpolate(const double* x, double hT) const;

private:
    std::string filename_for_time(double t) const;

    // internal helpers (templated on D=2/3, runtime dispatch from interpolate())
    template<unsigned int DD>
    InterpResult interpolate_impl(const double* x, double hT) const;
};

#endif
