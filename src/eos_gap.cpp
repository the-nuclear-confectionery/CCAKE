#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "../include/constants.h"
#include "../include/eos_gap.h"

using namespace constants;

////////////////////////////////////////////////////////////////////////////////
bool GapTable::load(const fs::path & path, bool normalize_by_T)
{
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "GapTable::load: cannot open " << path << std::endl;
        return false;
    }

    const double hc  = hbarc_MeVfm;
    const double hc1 = 1.0 / hc;
    const double hc2 = hc1 * hc1;
    const double hc3 = hc2 * hc1;
    const double hc4 = hc3 * hc1;

    // field indices in the row after the 4 axis columns:
    // p=0, s=1, B=2, S=3, Q=4, e=5, cs2=6, chiBB=7, ..., chiTT=16
    static_assert(N_GAP_THERMO == 17, "column count mismatch");

    int n_loaded = 0;
    int n_skipped_unstable = 0;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);
        double vals[4 + N_GAP_THERMO];
        bool ok = true;
        for (int i = 0; i < 4 + N_GAP_THERMO; ++i) {
            if (!(ss >> vals[i])) { ok = false; break; }
        }
        if (!ok) continue;

        GapPoint gp;
        // Apply axis scaling: MeV -> fm^-1
        double T_fm   = vals[0] / hc;
        double muB_fm = vals[1] / hc;
        double muQ_fm = vals[2] / hc;
        double muS_fm = vals[3] / hc;

        gp.T   = T_fm;
        gp.muB = muB_fm;
        gp.muQ = muQ_fm;
        gp.muS = muS_fm;

        // Apply field scaling (same as eos_table.cpp)
        if (normalize_by_T) {
            const double T2 = T_fm * T_fm;
            const double T3 = T2  * T_fm;
            const double T4 = T3  * T_fm;
            gp.thermo[0]  = vals[4]  * T4;  // p
            gp.thermo[1]  = vals[5]  * T3;  // s
            gp.thermo[2]  = vals[6]  * T3;  // B
            gp.thermo[3]  = vals[7]  * T3;  // S
            gp.thermo[4]  = vals[8]  * T3;  // Q
            gp.thermo[5]  = vals[9]  * T4;  // e
            gp.thermo[6]  = vals[10];        // cs2 (dimensionless)
            for (int k = 7; k < 17; ++k)
                gp.thermo[k] = vals[4 + k] * T2;  // chi's
        } else {
            gp.thermo[0]  = vals[4]  * hc4;  // p (MeV^4 -> fm^-4)
            gp.thermo[1]  = vals[5]  * hc3;  // s
            gp.thermo[2]  = vals[6]  * hc3;  // B
            gp.thermo[3]  = vals[7]  * hc3;  // S
            gp.thermo[4]  = vals[8]  * hc3;  // Q
            gp.thermo[5]  = vals[9]  * hc4;  // e
            gp.thermo[6]  = vals[10];          // cs2
            for (int k = 7; k < 17; ++k)
                gp.thermo[k] = vals[4 + k] * hc2;  // chi's
        }

        // Skip "truly unstable" rows: chiBB < 0 means the baryon susceptibility
        // is negative (genuine spinodal instability), which makes the derivative
        // matrix non-positive-definite and crashes downstream transport code.
        // Keep only the metastable subset of the gap (chiBB >= 0).
        //if (gp.thermo[7] < 0.0) {
        //    ++n_skipped_unstable;
        //    continue;
        //}

        // Mirror eos_table.cpp: clamp cs2 to >=0 (avoid negative-cs2 propagating
        // into transport/viscosity code that doesn't handle the unstable branch).
        //if (gp.thermo[6] < 0.0) gp.thermo[6] = 0.0;

        // Convenience lookup keys extracted from thermo
        gp.s    = gp.thermo[1];
        gp.rhoB = gp.thermo[2];
        gp.rhoS = gp.thermo[3];
        gp.rhoQ = gp.thermo[4];
        gp.e    = gp.thermo[5];

        points.push_back(gp);
        ++n_loaded;
    }

    loaded = (n_loaded > 0);
    std::cout << "GapTable::load: loaded " << n_loaded
              << " gap points from " << path
              << " (skipped " << n_skipped_unstable
              << " rows with chiBB<0)" << std::endl;
    return loaded;
}


////////////////////////////////////////////////////////////////////////////////
bool GapTable::lookup(double key_in, double rhoB_in, double rhoS_in, double rhoQ_in,
                      const std::string & mode, double tol,
                      double & T_out, double & muB_out, double & muQ_out, double & muS_out,
                      double thermo_out[N_GAP_THERMO],
                      bool use_energy) const
{
    if (points.empty()) return false;

    // Find two nearest points in (key, rhoB, rhoS, rhoQ) space
    double d1 = 1e300, d2 = 1e300;
    int    i1 = -1,    i2 = -1;

    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        const GapPoint & gp = points[i];
        const double key_gp = use_energy ? gp.e : gp.s;
        const double dk  = key_in  - key_gp;
        const double dB  = rhoB_in - gp.rhoB;
        const double dS  = rhoS_in - gp.rhoS;
        const double dQ  = rhoQ_in - gp.rhoQ;
        const double d   = std::sqrt(dk*dk + dB*dB + dS*dS + dQ*dQ);
        if (d < d1) { d2 = d1; i2 = i1; d1 = d; i1 = i; }
        else if (d < d2) { d2 = d; i2 = i; }
    }

    if (d1 > tol) return false;

    if (mode == "nearest" || i2 < 0) {
        // Return the single nearest point
        const GapPoint & gp = points[i1];
        T_out   = gp.T;
        muB_out = gp.muB;
        muQ_out = gp.muQ;
        muS_out = gp.muS;
        for (int k = 0; k < N_GAP_THERMO; ++k) thermo_out[k] = gp.thermo[k];
        return true;
    }

    // Linear interpolation: project target onto segment p1 -> p2
    const GapPoint & p1 = points[i1];
    const GapPoint & p2 = points[i2];

    const double key1 = use_energy ? p1.e : p1.s;
    const double key2 = use_energy ? p2.e : p2.s;

    // difference vector p2 - p1
    const double dk12 = key2 - key1;
    const double dB12 = p2.rhoB - p1.rhoB;
    const double dS12 = p2.rhoS - p1.rhoS;
    const double dQ12 = p2.rhoQ - p1.rhoQ;
    const double dot_dd = dk12*dk12 + dB12*dB12 + dS12*dS12 + dQ12*dQ12;

    double t = 0.0;
    if (dot_dd > 0.0) {
        // dot product of (target - p1) with (p2 - p1)
        const double dot_td = (key_in - key1)*dk12 + (rhoB_in - p1.rhoB)*dB12
                             + (rhoS_in - p1.rhoS)*dS12 + (rhoQ_in - p1.rhoQ)*dQ12;
        t = dot_td / dot_dd;
        t = std::max(0.0, std::min(1.0, t));
    }

    const double u = 1.0 - t;
    T_out   = u * p1.T   + t * p2.T;
    muB_out = u * p1.muB + t * p2.muB;
    muQ_out = u * p1.muQ + t * p2.muQ;
    muS_out = u * p1.muS + t * p2.muS;
    for (int k = 0; k < N_GAP_THERMO; ++k)
        thermo_out[k] = u * p1.thermo[k] + t * p2.thermo[k];

    return true;
}


////////////////////////////////////////////////////////////////////////////////
// Analytic override stub — returns false until filled in by the user.
bool GapTable::analytic_lookup(double /*key_in*/, double /*rhoB_in*/,
                               double /*rhoS_in*/, double /*rhoQ_in*/,
                               double & /*T_out*/, double & /*muB_out*/,
                               double & /*muQ_out*/, double & /*muS_out*/,
                               double /*thermo_out*/[N_GAP_THERMO]) const
{
    return false;
}
