#include "system_state_int.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>

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

    std::ostringstream fn;
    fn << output_dir << "/system_state_" << idx << ".dat";
    return fn.str();
}

void SystemState::load(double t)
{
    particles.clear();

    const std::string filename = filename_for_time(t);
    std::ifstream in(filename.c_str());
    if (!in)
        throw std::runtime_error("Cannot open " + filename);

    // first line: global time stored in file
    in >> file_time;

    while (true) {
        Particle p;

        if (!(in >> p.id >> p.t))
            break;

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

        // ---- SPH mass/weight ----
        // Your current file format does NOT show it.
        // So default to 1.0 for now (you can change this rule later).
        p.sph_mass = 1.0;

        // enforce dimension consistency
        if (D == 2) {
            p.r[2] = 0.0;
            p.u[2] = 0.0;
        }

        particles.push_back(p);
    }
}

// -------------------- interpolation --------------------

SystemState::InterpResult SystemState::interpolate(const double* x, double hT) const
{
    if (hT <= 0.0)
        throw std::runtime_error("SystemState::interpolate: hT must be > 0");

    if (D == 2) return interpolate_impl<2>(x, hT);
    if (D == 3) return interpolate_impl<3>(x, hT);

    throw std::runtime_error("SystemState::interpolate: D must be 2 or 3");
}

template<unsigned int DD>
SystemState::InterpResult SystemState::interpolate_impl(const double* x, double hT) const
{
    InterpResult out;
    //calculate sigmastar
    for (const auto& p : particles) {
        const double r = ccake::SPHkernel<DD>::distance(p.r, x);
        if (r >= 2.0 * hT) continue;

        const double W = ccake::SPHkernel<DD>::kernel(r, hT);
        const double m = p.sph_mass;
        out.sigma_star += m * W;
    }


    for (const auto& p : particles) {
        const double r = ccake::SPHkernel<DD>::distance(p.r, x);
        if (r >= 2.0 * hT) continue;

        const double W = ccake::SPHkernel<DD>::kernel(r, hT);
        const double m = p.sph_mass;
        out.e     += m * p.e    * W/out.sigma_star;
        out.s     += m * p.s    * W/out.sigma_star;
        out.rhoB  += m * p.rhoB * W/out.sigma_star;
        out.rhoS  += m * p.rhoS * W/out.sigma_star;
        out.rhoQ  += m * p.rhoQ * W/out.sigma_star;
        out.T     += m * p.T    * W/out.sigma_star;
    }

    return out;
}
