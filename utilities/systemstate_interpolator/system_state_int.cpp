#include "system_state_int.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <regex>

#include <Kokkos_Core.hpp>

namespace fs = std::filesystem;

// ---- internal helpers -------------------------------------------------------

namespace {

SystemState::InterpResult buildResult(
    double mW, double meW, double msW,
    double mrBW, double mrSW, double mrQW, double mTW)
{
    SystemState::InterpResult out;
    out.sigma_star = mW;
    if (mW > 0.0) {
        out.e    = meW  / mW;
        out.s    = msW  / mW;
        out.rhoB = mrBW / mW;
        out.rhoS = mrSW / mW;
        out.rhoQ = mrQW / mW;
        out.T    = mTW  / mW;
    }
    return out;
}

SystemState::InterpResult blend(
    const SystemState::InterpResult& a,
    const SystemState::InterpResult& b,
    double alpha)
{
    const double w = 1.0 - alpha;
    SystemState::InterpResult r;
    r.sigma_star = w * a.sigma_star + alpha * b.sigma_star;
    r.e          = w * a.e          + alpha * b.e;
    r.s          = w * a.s          + alpha * b.s;
    r.rhoB       = w * a.rhoB       + alpha * b.rhoB;
    r.rhoS       = w * a.rhoS       + alpha * b.rhoS;
    r.rhoQ       = w * a.rhoQ       + alpha * b.rhoQ;
    r.T          = w * a.T          + alpha * b.T;
    return r;
}

SystemState::SphFields sphFrom(const SystemState::Particle& p)
{
    SystemState::SphFields f;
    f.r[0] = p.r[0]; f.r[1] = p.r[1]; f.r[2] = p.r[2];
    f.e    = p.e;    f.s    = p.s;
    f.rhoB = p.rhoB; f.rhoS = p.rhoS; f.rhoQ = p.rhoQ;
    f.T    = p.T;    f.sph_mass = p.sph_mass;
    return f;
}

} // namespace

// -------------------- constructor --------------------

SystemState::SystemState(int D_from_config,
                         const std::string& output_dir_,
                         double t0_, double dt_,
                         int output_steps_)
    : D(D_from_config),
      t0(t0_),
      dt(dt_),
      output_steps(output_steps_),
      output_dir(output_dir_),
      file_time(0.0)
{
    if (D != 2 && D != 3)
        throw std::runtime_error("SystemState: D must be 2 or 3");
    if (dt <= 0.0)
        throw std::runtime_error("SystemState: dt must be > 0");
    if (output_steps <= 0)
        throw std::runtime_error("SystemState: output_steps must be > 0");
}

// -------------------- filename helpers --------------------

std::string SystemState::filename_for_index(long long idx) const
{
    return output_dir + "/system_state_" + std::to_string(idx) + ".dat";
}

std::string SystemState::filename_for_time(double t) const
{
    const double stride = dt * double(output_steps);
    const double x = (t - t0) / stride;
    const long long idx = llround(x);

    const double t_expected = t0 + idx * stride;
    const double tol = 1e-10 * (1.0 + std::fabs(t_expected));
    if (std::fabs(t - t_expected) > tol) {
        std::ostringstream oss;
        oss << "Requested t=" << t
            << " not on output grid (nearest "
            << t_expected << ")";
        throw std::runtime_error(oss.str());
    }

    if (idx < 0)
        throw std::runtime_error("Requested t < t0");

    return filename_for_index(idx);
}

// -------------------- internal file loader --------------------

void SystemState::load_into(const std::string& filename,
                             TimeStep& ts, double h_cell) const
{
    std::ifstream in(filename.c_str());
    if (!in)
        throw std::runtime_error("Cannot open " + filename);

    in >> ts.time;
    ts.particles.clear();

    while (true) {
        Particle p;
        if (!(in >> p.id >> p.t)) break;

        in >> p.r[0] >> p.r[1] >> p.r[2];
        in >> p.p >> p.T >> p.muB >> p.muS >> p.muQ;
        in >> p.e >> p.rhoB >> p.rhoS >> p.rhoQ >> p.s;
        in >> p.eta_pi >> p.zeta_Pi >> p.tau_Pi >> p.tau_pi;
        in >> p.theta;
        in >> p.invRe_shear >> p.invRe_bulk;
        in >> p.kn_shear >> p.kn_bulk;

        for (int i = 0; i < 9; ++i)
            in >> p.shv[i];

        in >> p.u[0] >> p.u[1] >> p.u[2];
        in >> p.gamma;
        in >> p.Freeze;
        in >> p.eos;
        in >> p.causality;

        p.sph_mass = 1.0;  // not in file; all particles carry equal SPH mass

        if (D == 2) { p.r[2] = 0.0; p.u[2] = 0.0; }

        ts.particles.push_back(p);
    }

    if (ts.particles.empty()) return;

    // ---- build Cabana LinkedCellList ----
    const int N = (int)ts.particles.size();
    ts.cell_size = h_cell;

    // pos_view is local: only needed during CellList construction
    Kokkos::View<double*[3], Kokkos::HostSpace> pos_view(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "positions"), N);

    // fill positions and compute bounding box in one pass
    double lo[3] = { 1e30,  1e30,  1e30};
    double hi[3] = {-1e30, -1e30, -1e30};
    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) {
            pos_view(i, d) = ts.particles[i].r[d];
            lo[d] = std::min(lo[d], ts.particles[i].r[d]);
            hi[d] = std::max(hi[d], ts.particles[i].r[d]);
        }
    }

    // For 2D: collapse z to a single cell so all particles land in cell k=0
    if (D == 2) {
        lo[2] = -0.5 * h_cell;
        hi[2] =  0.5 * h_cell;
    }

    // Add one cell of padding so boundary particles are included
    for (int d = 0; d < 3; ++d) {
        lo[d] -= h_cell;
        hi[d] += h_cell;
        ts.lo[d] = lo[d];
    }

    double delta[3] = {h_cell, h_cell, h_cell};
    ts.cell_list = CellList(pos_view, delta, lo, hi);
    // pos_view destroyed here; CellList owns all needed data

    // ---- build cell-sorted compact SPH array (Opt 2+3) ----
    // Store SphFields in permuted order: ts.sph[s] corresponds directly to
    // cell bin offset s, eliminating permutation() indirection in the hot loop.
    ts.sph.resize(N);
    for (int s = 0; s < N; ++s)
        ts.sph[s] = sphFrom(ts.particles[(int)ts.cell_list.permutation(s)]);
}

// -------------------- load (legacy single-snapshot) --------------------

void SystemState::load(double t)
{
    particles.clear();
    TimeStep ts;
    load_into(filename_for_time(t), ts, 1.0);
    file_time = ts.time;
    particles  = std::move(ts.particles);
}

// -------------------- load_all --------------------

void SystemState::load_all(double h_cell)
{
    if (!fs::is_directory(output_dir))
        throw std::runtime_error("load_all: not a directory: " + output_dir);

    static const std::regex pat("system_state_(\\d+)\\.dat");
    std::vector<std::pair<long long, fs::path>> entries;
    for (const auto& entry : fs::directory_iterator(output_dir)) {
        std::smatch m;
        std::string name = entry.path().filename().string();
        if (std::regex_match(name, m, pat))
            entries.emplace_back(std::stoll(m[1].str()), entry.path());
    }

    if (entries.empty())
        throw std::runtime_error(
            "load_all: no system_state_*.dat files in " + output_dir);

    std::sort(entries.begin(), entries.end());

    evolution.clear();
    evolution.reserve(entries.size());
    for (const auto& [idx, path] : entries) {
        TimeStep ts;
        load_into(path.string(), ts, h_cell);
        evolution.push_back(std::move(ts));
    }
    last_lo_idx_ = 0;
}

// -------------------- SPH interpolation (cell list) --------------------

template<unsigned int DD>
SystemState::InterpResult SystemState::interpolate_impl(
    const TimeStep& ts, const double* x, double hT) const
{
    const double R  = 2.0 * hT;
    const double R2 = R * R;  // precompute for r² early-exit (Opt 1)
    const auto& lcl = ts.cell_list;

    // cell index of query point
    int qi[3];
    for (int d = 0; d < 3; ++d) {
        qi[d] = (int)std::floor((x[d] - ts.lo[d]) / ts.cell_size);
        qi[d] = std::max(0, std::min(lcl.numBin(d) - 1, qi[d]));
    }

    const int cr = (int)std::ceil(R / ts.cell_size);

    const int i0 = std::max(0, qi[0]-cr), i1 = std::min(lcl.numBin(0)-1, qi[0]+cr);
    const int j0 = std::max(0, qi[1]-cr), j1 = std::min(lcl.numBin(1)-1, qi[1]+cr);
    const int k0 = (DD==2) ? 0 : std::max(0, qi[2]-cr);
    const int k1 = (DD==2) ? 0 : std::min(lcl.numBin(2)-1, qi[2]+cr);

    double mW=0, meW=0, msW=0, mrBW=0, mrSW=0, mrQW=0, mTW=0;

    for (int ci = i0; ci <= i1; ++ci)
    for (int cj = j0; cj <= j1; ++cj)
    for (int ck = k0; ck <= k1; ++ck) {
        const auto off = lcl.binOffset(ci, cj, ck);
        const auto cnt = lcl.binSize(ci, cj, ck);
        // ts.sph[] is cell-sorted: ts.sph[off..off+cnt-1] are the particles in
        // this bin — no permutation() call needed (Opt 3)
        for (auto s = off; s < off + (decltype(off))cnt; ++s) {
            const SphFields& f = ts.sph[s];
            // r² early-exit: skip sqrt for rejected particles (Opt 1)
            double r2 = 0.0;
            for (unsigned d = 0; d < DD; ++d) {
                const double diff = f.r[d] - x[d];
                r2 += diff * diff;
            }
            if (r2 >= R2) continue;
            const double r = std::sqrt(r2);
            const double w = f.sph_mass * ccake::SPHkernel<DD>::kernel(r, hT);
            mW   += w;
            meW  += f.e    * w;
            msW  += f.s    * w;
            mrBW += f.rhoB * w;
            mrSW += f.rhoS * w;
            mrQW += f.rhoQ * w;
            mTW  += f.T    * w;
        }
    }

    return buildResult(mW, meW, msW, mrBW, mrSW, mrQW, mTW);
}

// -------------------- SPH interpolation (legacy, no cell list) --------------------

template<unsigned int DD>
SystemState::InterpResult SystemState::interpolate_impl_slow(
    const std::vector<Particle>& parts,
    const double* x, double hT) const
{
    const double R2 = 4.0 * hT * hT;  // (2*hT)² — r² early-exit (Opt 1)
    double mW=0, meW=0, msW=0, mrBW=0, mrSW=0, mrQW=0, mTW=0;

    for (const auto& p : parts) {
        double r2 = 0.0;
        for (unsigned d = 0; d < DD; ++d) {
            const double diff = p.r[d] - x[d];
            r2 += diff * diff;
        }
        if (r2 >= R2) continue;
        const double r = std::sqrt(r2);
        const double w = p.sph_mass * ccake::SPHkernel<DD>::kernel(r, hT);
        mW   += w;
        meW  += p.e    * w;
        msW  += p.s    * w;
        mrBW += p.rhoB * w;
        mrSW += p.rhoS * w;
        mrQW += p.rhoQ * w;
        mTW  += p.T    * w;
    }

    return buildResult(mW, meW, msW, mrBW, mrSW, mrQW, mTW);
}

// public dispatch (legacy, uses this->particles)
SystemState::InterpResult SystemState::interpolate(
    const double* x, double hT) const
{
    if (hT <= 0.0)
        throw std::runtime_error("interpolate: hT must be > 0");
    if (D == 2) return interpolate_impl_slow<2>(particles, x, hT);
    if (D == 3) return interpolate_impl_slow<3>(particles, x, hT);
    throw std::runtime_error("interpolate: D must be 2 or 3");
}

// -------------------- interpolate_at (time + space) --------------------

// Core implementation given a known bracket index — no mutable state, thread-safe.
SystemState::InterpResult SystemState::interpolate_at_impl(
    double t, const double* x, double hT, std::size_t lo_idx) const
{
    const TimeStep& low  = evolution[lo_idx];
    const TimeStep& high = evolution[lo_idx + 1];

    const double tol = 1e-12 * (1.0 + std::fabs(low.time));
    if (std::fabs(t - low.time) < tol)
        return D==2 ? interpolate_impl<2>(low, x, hT) : interpolate_impl<3>(low, x, hT);

    const double alpha = (t - low.time) / (high.time - low.time);

    if (D == 2)
        return blend(interpolate_impl<2>(low, x, hT), interpolate_impl<2>(high, x, hT), alpha);
    return blend(interpolate_impl<3>(low, x, hT), interpolate_impl<3>(high, x, hT), alpha);
}

// Helper: binary search for the lower bracket index.
static std::size_t findBracket(const std::vector<SystemState::TimeStep>& ev, double t)
{
    auto it = std::upper_bound(ev.begin(), ev.end(), t,
        [](double val, const SystemState::TimeStep& ts){ return val < ts.time; });
    // clamp to valid interior range [0, N-2]
    if (it == ev.end())    return ev.size() - 2;
    if (it == ev.begin())  return 0;
    return (std::size_t)std::distance(ev.begin(), std::prev(it));
}

// Sequential version: uses bracket cache (Opt 5). Not thread-safe.
SystemState::InterpResult SystemState::interpolate_at(
    double t, const double* x, double hT) const
{
    if (evolution.empty())
        throw std::runtime_error("interpolate_at: call load_all() first");
    if (hT <= 0.0)
        throw std::runtime_error("interpolate_at: hT must be > 0");

    const double t_min = evolution.front().time;
    const double t_max = evolution.back().time;
    if (t < t_min || t > t_max) {
        std::ostringstream oss;
        oss << "interpolate_at: t=" << t
            << " out of range [" << t_min << ", " << t_max << "]";
        throw std::runtime_error(oss.str());
    }

    // edge cases: t exactly at first or last snapshot
    if (t <= t_min)
        return D==2 ? interpolate_impl<2>(evolution.front(), x, hT)
                    : interpolate_impl<3>(evolution.front(), x, hT);
    if (t >= t_max)
        return D==2 ? interpolate_impl<2>(evolution.back(), x, hT)
                    : interpolate_impl<3>(evolution.back(), x, hT);

    // Bracket cache (Opt 5): check if t is in the last-used bracket first
    const std::size_t cached = last_lo_idx_;
    if (cached + 1 < evolution.size() &&
        t >= evolution[cached].time && t <= evolution[cached+1].time)
    {
        return interpolate_at_impl(t, x, hT, cached);
    }

    // Cache miss: binary search and update cache
    const std::size_t lo_idx = findBracket(evolution, t);
    last_lo_idx_ = lo_idx;
    return interpolate_at_impl(t, x, hT, lo_idx);
}

// Parallel batch version (Opt 4): thread-safe, no shared mutable state.
std::vector<SystemState::InterpResult> SystemState::interpolate_many(
    const std::vector<Query>& queries) const
{
    if (evolution.empty())
        throw std::runtime_error("interpolate_many: call load_all() first");

    const int M = (int)queries.size();
    std::vector<InterpResult> out(M);

    Kokkos::parallel_for("SystemState::interpolate_many",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, M),
        [&](int i) {
            const Query& q = queries[i];
            const double t_min = evolution.front().time;
            const double t_max = evolution.back().time;

            if (q.t <= t_min) {
                out[i] = D==2 ? interpolate_impl<2>(evolution.front(), q.x, q.hT)
                               : interpolate_impl<3>(evolution.front(), q.x, q.hT);
                return;
            }
            if (q.t >= t_max) {
                out[i] = D==2 ? interpolate_impl<2>(evolution.back(), q.x, q.hT)
                               : interpolate_impl<3>(evolution.back(), q.x, q.hT);
                return;
            }

            const std::size_t lo = findBracket(evolution, q.t);
            out[i] = interpolate_at_impl(q.t, q.x, q.hT, lo);
        });
    Kokkos::fence();
    return out;
}
