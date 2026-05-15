#ifndef EOS_GAP_H
#define EOS_GAP_H

#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static constexpr int N_GAP_THERMO = 17;

/// A single scattered point in the auxiliary gap table.
/// Lookup keys are the conserved densities (fm^-3, after unit scaling).
/// The 17 thermo fields follow the same ordering as the main EoS table:
/// p, s, B, S, Q, e, cs2, chiBB, chiQQ, chiSS, chiBQ, chiBS, chiQS, chiTB, chiTQ, chiTS, chiTT
struct GapPoint {
    double s, rhoB, rhoS, rhoQ; ///< lookup keys in entropy/charge-density space (fm^-3)
    double e;                   ///< energy density (fm^-4) — alternate lookup key
    double T, muB, muQ, muS;   ///< phase-diagram coordinates (fm^-1, after hbarc scaling)
    double thermo[N_GAP_THERMO];///< full thermodynamics array (after unit scaling)
};

/// Auxiliary scattered-points table for the gap region of a first-order transition.
/// Stores rows dropped by the stable-phase deduplication of the main EoS table so
/// that the rootfinder can fall back to them when a target (s,rhoB,...) lies inside
/// the spinodal gap.
class GapTable {
public:
    std::vector<GapPoint> points;
    bool loaded = false;

    /// Read scattered-points gap file (same column format as main table, no grid-size header).
    /// Applies the same unit conversions as EoS_table uses for the main table.
    bool load(const fs::path & path, bool normalize_by_T);

    /// Find the nearest gap point(s) and return interpolated (T,mu,thermo).
    /// @param use_energy  if true, use energy density as first lookup key instead of entropy density
    /// @return true if a point within tolerance was found
    bool lookup(double key_in, double rhoB_in, double rhoS_in, double rhoQ_in,
                const std::string & mode, double tol,
                double & T_out, double & muB_out, double & muQ_out, double & muS_out,
                double thermo_out[N_GAP_THERMO],
                bool use_energy = false) const;

    /// Hard-coded analytic override stub.
    /// Returns false (stub) — fill in with EoS-specific T_c(muB) fits when available.
    bool analytic_lookup(double key_in, double rhoB_in, double rhoS_in, double rhoQ_in,
                         double & T_out, double & muB_out, double & muQ_out, double & muS_out,
                         double thermo_out[N_GAP_THERMO]) const;
};

#endif
