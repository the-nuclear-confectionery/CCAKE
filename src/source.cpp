#include <iomanip>

#include "source.h"

namespace ccake
{
// Source<D> is held by SPHWorkstation<D>, instantiated for D = 1, 2, 3.
template class Source<1>;
template class Source<2>;
template class Source<3>;

namespace
{
// Charges of an AMPT parton. AMPT writes PDG flavor codes (1=d, 2=u, 3=s, ...,
// negative for antiquarks; checked against its masses and Au's Z/A).
void ampt_parton_charges(int pid, double& baryon_charge, double& strangeness_charge,
                         double& electric_charge)
{
  baryon_charge = 0.0;
  strangeness_charge = 0.0;
  electric_charge = 0.0;
  if (pid == 21) return;   // gluon

  int flavor = pid;
  double sign = 1.0;       // antiquarks flip every charge
  if (pid < 0)
  {
    flavor = -pid;
    sign = -1.0;
  }
  if (flavor < 1 || flavor > 6) throw std::runtime_error("Unknown PID in source file");

  baryon_charge = sign/3.0;
  if (flavor == 2 || flavor == 4 || flavor == 6)   // u, c, t
    electric_charge = sign*(2.0/3.0);
  else                                             // d, s, b
    electric_charge = sign*(-1.0/3.0);
  if (flavor == 3) strangeness_charge = -sign;
}

// The whitespace-separated numbers on one line of a source file.
std::vector<double> numbers_in(const std::string& line)
{
  std::vector<double> numbers;
  std::istringstream stream(line);
  double value;
  while (stream >> value) numbers.push_back(value);
  return numbers;
}
}

//==============================================================================
template <unsigned int D>
void Source<D>::initialize( std::shared_ptr<Settings> settingsPtr_in,
                            std::shared_ptr<SystemState<D>> systemPtr_in )
{
  settingsPtr = settingsPtr_in;
  systemPtr   = systemPtr_in;
  formatted_output::report("Initializing source terms");

  const std::string model = settingsPtr->source_model;
  if (model == "BBMG")
  {
    has_bbmg_sources = true;
    formatted_output::detail("Initializing real-time BBMG jet source term");
    return;
  }
  if (model != "AMPT" && model != "SMASH" && model != "MCGLB")
  {
    std::cerr << "ERROR: unknown hydro.source.model '" << model
              << "' (AMPT, SMASH, MCGLB or BBMG)." << std::endl;
    abort();
  }
  if (settingsPtr->source_type != "ic")
  {
    std::cerr << "ERROR: hydro.source.model " << model << " needs hydro.source.type: ic." << std::endl;
    abort();
  }

  // Read once; advance_sources() then releases each budget step by step.
  std::vector<SourceRecord> records;
  if (model == "AMPT")  read_ampt_file(records);
  if (model == "SMASH") read_smash_file(records);
  if (model == "MCGLB") read_mcglb_file(records);
  load_sources(records);
  has_ic_sources = true;

  std::string deposit_mode = "instant deposit";
  if (settingsPtr->source_propagate) deposit_mode = "propagate & deplete (constant-rate loss)";
  formatted_output::detail("Loaded " + std::to_string(ic_sources.size()) + " " + model
                           + " sources: " + deposit_mode);
}

//==============================================================================
template <unsigned int D>
void Source<D>::deposit_sources()
{
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  double t  = systemPtr->t;
  double dt = settingsPtr->dt;

  if (has_bbmg_sources) build_bbmg_sources(t);
  if (has_ic_sources)   build_step_sources();

  double smearing_radius     = settingsPtr->smearing_radius;      // transverse
  double smearing_radius_eta = settingsPtr->smearing_radius_eta;  // longitudinal (eta)

  // reset source terms
  auto reset_source_terms = KOKKOS_LAMBDA(const int is, const int ia)
  {
    device_hydro_scalar.access(is, ia, ccake::hydro_info::j0_ext) = 0.0;
    for (int idir = 0; idir < D; ++idir)
      device_hydro_vector.access(is, ia, ccake::hydro_info::j_ext, idir) = 0.0;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::rho_B_ext) = 0.0;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::rho_S_ext) = 0.0;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::rho_Q_ext) = 0.0;
  };
  Cabana::simd_parallel_for(simd_policy, reset_source_terms, "reset_source_terms");
  Kokkos::fence();

  // Views of the source AoSoA
  SRC_VIEW(s_, sources);

  // One kernel for all sources: each SPH particle loops over them.
  const double normalization        = settingsPtr->source_normalization;
  const double inv_ntest            = 1.0/settingsPtr->source_test_particles; // SMASH test-particle norm
  const bool   is_hyperbolic         = (settingsPtr->coordinate_system == "hyperbolic");
  const bool   deposit_baryon        = settingsPtr->baryon_source;
  const bool   deposit_strangeness   = settingsPtr->strangeness_source;
  const bool   deposit_electric      = settingsPtr->electric_source;
  const int    num_sources           = (int)n_sources;

  // host-side count of sources depositing this step (for the log)
  int num_sources_deposited = 0;
  for (int isrc = 0; isrc < num_sources; ++isrc)
    if (s_energy(isrc) > 0.0) num_sources_deposited++;

  auto deposit_all = KOKKOS_LAMBDA(const int is, const int ia)
  {
    double particle_position[D];
    for (int idir = 0; idir < D; ++idir)
      particle_position[idir] = device_position.access(is, ia, idir);

    double energy_current = 0.0, momentum_current[D];
    double rho_baryon = 0.0, rho_strangeness = 0.0, rho_electric = 0.0;
    for (int idir = 0; idir < D; ++idir) momentum_current[idir] = 0.0;

    for (int isrc = 0; isrc < num_sources; ++isrc)
    {
      if (s_energy(isrc) <= 0.0) continue;                 // nothing to deposit this step

      double source_position[D];
      for (int idir = 0; idir < D; ++idir) source_position[idir] = s_position(isrc, idir);
      // factorized kernel: transverse smearing_radius x longitudinal smearing_radius_eta
      double kernel_value = SPHkernel<D>::kernel(particle_position, source_position,
                                                 smearing_radius, smearing_radius_eta);
      if (kernel_value <= 0.0) continue;
      kernel_value /= s_renormalization(isrc);             // Shepard: 1/coverage, capped (1 = none)
      if (is_hyperbolic) kernel_value /= t;                // Milne volume Jacobian

      const double weight = inv_ntest * normalization * kernel_value / dt; // -> Energy/time, Momentum/time
      energy_current += weight * s_energy(isrc);
      for (int idir = 0; idir < D; ++idir) momentum_current[idir] += weight * s_momentum(isrc, idir);
      if (deposit_baryon)      rho_baryon      += inv_ntest * s_baryon_charge(isrc)      * kernel_value / dt;
      if (deposit_strangeness) rho_strangeness += inv_ntest * s_strangeness_charge(isrc) * kernel_value / dt;
      if (deposit_electric)    rho_electric    += inv_ntest * s_electric_charge(isrc)    * kernel_value / dt;
    }

    device_hydro_scalar.access(is, ia, ccake::hydro_info::j0_ext) += energy_current;
    for (int idir = 0; idir < D; ++idir)
      device_hydro_vector.access(is, ia, ccake::hydro_info::j_ext, idir) += momentum_current[idir];
    device_hydro_scalar.access(is, ia, ccake::hydro_info::rho_B_ext) += rho_baryon;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::rho_S_ext) += rho_strangeness;
    device_hydro_scalar.access(is, ia, ccake::hydro_info::rho_Q_ext) += rho_electric;
  };
  Cabana::simd_parallel_for(simd_policy, deposit_all, "deposit_all_sources");
  Kokkos::fence();

  formatted_output::detail("Added source terms for " + std::to_string(num_sources_deposited)
                           + " source particles out of " + std::to_string(n_sources));
}

//==============================================================================
template <unsigned int D>
void Source<D>::transform_to_hyperbolic( bool is_hyperbolic,
    double& px, double& py, double& pz,
    double& x,  double& y,  double& z,
    double& t,  double& E )
{
  if (!is_hyperbolic)
    return;

  if (t*t < z*z)
    return;

  double tau   = std::sqrt(t*t - z*z);
  double eta_s = 0.5 * std::log((t + z) / (t - z));

  double p_tau =  E * std::cosh(eta_s) - pz * std::sinh(eta_s);
  double p_eta = (-E * std::sinh(eta_s) + pz * std::cosh(eta_s))/tau;

  // Overwrite momentum & coordinates
  pz = p_eta;
  E  = p_tau;
  t  = tau;
  z  = eta_s;
}

//==============================================================================
template <unsigned int D>
bool Source<D>::passes_time_cut(double t, double z) const
{
  double start_time = t;
  if (settingsPtr->coordinate_system == "hyperbolic")
  {
    // only keep particles in the light cone
    if (t*t < z*z)
    {
      formatted_output::detail("Skipping source entry with t^2 < z^2 at t = "
                               + std::to_string(t) + " and z = " + std::to_string(z));
      return false;
    }
    start_time = std::sqrt(t*t - z*z); // use tau for the comparison in hyperbolic coordinates
  }
  return start_time > settingsPtr->t0;
}

//==============================================================================
template <unsigned int D>
void Source<D>::load_sources(const std::vector<SourceRecord>& records)
{
  const bool is_hyperbolic = (settingsPtr->coordinate_system == "hyperbolic");

  ic_sources.clear();
  for (const SourceRecord& record : records)
  {
    double px = record.px, py = record.py, pz = record.pz, E = record.E;
    double x = record.x, y = record.y, z = record.z, t = record.t;

    // deposit_at_rapidity: place the source at eta_s = its momentum rapidity, at the
    // same formation proper time (low-energy partons form near z = 0 with |y| up to 2).
    if (settingsPtr->source_deposit_at_rapidity && is_hyperbolic)
    {
      double energy = E;
      if (energy <= 0.0) energy = std::sqrt(px*px + py*py + pz*pz + record.mass*record.mass);
      if (energy > std::abs(pz))
      {
        const double rapidity = 0.5*std::log((energy + pz)/(energy - pz));
        const double tau_form = std::sqrt(std::max(t*t - z*z, 0.0));
        t = tau_form*std::cosh(rapidity);
        z = tau_form*std::sinh(rapidity);
      }
    }

    // Position to Milne; momentum stays Cartesian (the Milne velocity is rebuilt
    // each step). The budget is the Milne energy p^tau.
    transform_to_hyperbolic(is_hyperbolic, px, py, pz, x, y, z, t, E);

    ICSource source;
    source.px = record.px/hbarc_GeVfm;
    source.py = record.py/hbarc_GeVfm;
    source.pz = record.pz/hbarc_GeVfm;
    source.mass = record.mass/hbarc_GeVfm;
    source.E0   = E/hbarc_GeVfm;
    source.baryon  = record.baryon_charge;
    source.strange = record.strangeness_charge;
    source.charge  = record.electric_charge;
    source.formation_time = t;
    source.PID    = record.pid;
    source.type   = record.type;
    source.parent = record.parent;
    source.x = x;
    source.y = y;
    source.z = z;   // eta_s
    source.remaining = source.E0;
    ic_sources.push_back(source);
  }

  n_sources = ic_sources.size();
  sources   = source_array("sources", n_sources);
}

//==============================================================================
// AMPT parton file: one header line, then  pid px py pz mass x y z t ...
template <unsigned int D>
void Source<D>::read_ampt_file(std::vector<SourceRecord>& records) const
{
  std::ifstream infile(settingsPtr->source_input_file);
  if (!infile.is_open())
    throw std::runtime_error("Failed to open source file: " + settingsPtr->source_input_file);

  std::string line;
  std::getline(infile, line);   // header
  while (std::getline(infile, line))
  {
    const std::vector<double> column = numbers_in(line);
    if (column.size() < 12) continue;

    SourceRecord record;
    record.pid  = (int)column[0];
    record.px   = column[1];
    record.py   = column[2];
    record.pz   = column[3];
    record.mass = column[4];
    record.x    = column[5];
    record.y    = column[6];
    record.z    = column[7];
    record.t    = column[8];
    if (!passes_time_cut(record.t, record.z)) continue;
    record.E = std::sqrt(record.px*record.px + record.py*record.py
                       + record.pz*record.pz + record.mass*record.mass);
    ampt_parton_charges(record.pid, record.baryon_charge, record.strangeness_charge,
                        record.electric_charge);
    records.push_back(record);
  }
}

//==============================================================================
// SMASH OSCAR2013Extended particle list: t x y z mass p0 px py pz pdg id charge ...
// with baryon number and strangeness in columns 21 and 22. The enable_* flags decide
// which charges are deposited.
template <unsigned int D>
void Source<D>::read_smash_file(std::vector<SourceRecord>& records) const
{
  std::ifstream infile(settingsPtr->source_input_file);
  if (!infile.is_open())
    throw std::runtime_error("Failed to open source file: " + settingsPtr->source_input_file);

  std::string line;
  while (std::getline(infile, line))
  {
    if (line.empty() || line[0] == '#') continue;
    const std::vector<double> column = numbers_in(line);
    if (column.size() < 22) continue;

    SourceRecord record;
    record.t    = column[0];
    record.x    = column[1];
    record.y    = column[2];
    record.z    = column[3];
    record.mass = column[4];
    record.px   = column[6];
    record.py   = column[7];
    record.pz   = column[8];
    record.pid  = (int)column[9];
    record.electric_charge    = column[11];
    record.baryon_charge      = column[20];
    record.strangeness_charge = column[21];
    if (!passes_time_cut(record.t, record.z)) continue;
    record.E = std::sqrt(record.px*record.px + record.py*record.py
                       + record.pz*record.pz + record.mass*record.mass);
    records.push_back(record);
  }
}

//==============================================================================
// 3dMCGlauber string elements from utilities/mcglb2ccake.py --source-file:
//   t x y z E px py pz baryon strange electric [type parent]
// Cartesian formation point [fm], 4-momentum [GeV], explicit (possibly fractional)
// charges. Optional tags: type 0 body / 1 left / 2 right remnant / 3 carrier, parent row.
template <unsigned int D>
void Source<D>::read_mcglb_file(std::vector<SourceRecord>& records) const
{
  std::ifstream infile(settingsPtr->source_input_file);
  if (!infile.is_open())
    throw std::runtime_error("Failed to open source file: " + settingsPtr->source_input_file);

  std::string line;
  while (std::getline(infile, line))
  {
    if (line.empty() || line[0] == '#') continue;
    const std::vector<double> column = numbers_in(line);
    if (column.size() < 11) continue;

    SourceRecord record;
    record.t  = column[0];
    record.x  = column[1];
    record.y  = column[2];
    record.z  = column[3];
    record.E  = column[4];
    record.px = column[5];
    record.py = column[6];
    record.pz = column[7];
    record.baryon_charge      = column[8];
    record.strangeness_charge = column[9];
    record.electric_charge    = column[10];
    if (column.size() > 11) record.type   = (int)column[11];
    if (column.size() > 12) record.parent = (int)column[12];
    if (!passes_time_cut(record.t, record.z)) continue;

    // invariant mass of the element (its transverse mass; strings carry no pT)
    const double p_sq = record.px*record.px + record.py*record.py + record.pz*record.pz;
    const double m_sq = record.E*record.E - p_sq;
    record.mass = 0.0;
    if (m_sq > 0.0) record.mass = std::sqrt(m_sq);
    record.pid = 99;
    records.push_back(record);
  }
}

//==============================================================================
template <unsigned int D>
void Source<D>::build_bbmg_sources(double t)
{
  // Each jet's per-step energy loss, with its momentum along the jet direction
  // (jets are already in simulation coordinates).
  if (bbmgPtr == nullptr)
  {
    formatted_output::detail("BBMG source requested but no BBMG object was set.");
    return;
  }

  bbmgPtr->copy_device_to_host_BBMG();
  auto& jets = bbmgPtr->jetInfo_host;

  n_sources = jets.size();
  if (n_sources == 0) return;
  if (sources.size() != n_sources) sources.resize(n_sources);
  SRC_VIEW(s_, sources);

  double tau = t;
  for (std::size_t i = 0; i < jets.size(); ++i)
  {
    double dE = jets[i].e_loss;

    // unit Milne velocity, as in BBMG::propagate: the deposit is massless, E = |p|
    double delta_y = 0.0;
    if (D == 3) delta_y = jets[i].rapidity - jets[i].eta;
    double sech_delta_y = 1.0 / std::cosh(delta_y);
    double direction[3] = { std::cos(jets[i].phi) * sech_delta_y,
                            std::sin(jets[i].phi) * sech_delta_y,
                            0.0 };
    if (D == 3) direction[2] = std::tanh(delta_y) / tau;

    for (int idir = 0; idir < 3; ++idir)
    {
      s_position(i, idir) = 0.0;
      if (idir < (int)D) s_position(i, idir) = jets[i].r[idir];
      s_momentum(i, idir) = dE * direction[idir]; // massless: |p|_Milne = E along direction
    }
    s_energy(i) = dE;
    s_baryon_charge(i) = 0.0;
    s_strangeness_charge(i) = 0.0;
    s_electric_charge(i) = 0.0;
    s_renormalization(i) = 1.0;
  }
}

//==============================================================================
template <unsigned int D>
double Source<D>::delivered_fraction(double coverage) const
{
  const double max_renormalization_factor = std::max(1.0, settingsPtr->source_max_renormalization_factor);
  return std::min(1.0, coverage * max_renormalization_factor);
}

//==============================================================================
// kernel_renormalize: each active source's kernel coverage at its current position.
// coverage = sum_a W(r_a - r_src) m_a / sigma_a is the discrete kernel integral (1 where
// the fluid covers the kernel); deposits are divided by it, and sources below min_coverage
// wait, budget intact, until the fluid reaches them. The same sum over frozen particles
// records the deposits that never reach Cooper-Frye (record_frozen_deposit).
template <unsigned int D>
void Source<D>::compute_coverage()
{
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  const int    num_particles       = (int)systemPtr->cabana_particles.size();
  const double t                   = systemPtr->t;
  const double smearing_radius     = settingsPtr->smearing_radius;
  const double smearing_radius_eta = settingsPtr->smearing_radius_eta;

  for (ICSource& source : ic_sources)
  {
    source.coverage = 0.0;
    source.frozen_weight = 0.0;
    source.permanently_frozen_weight = 0.0;
    source.frozen_at_birth_particles.clear();
    source.frozen_at_birth_weights.clear();
    if (t < source.formation_time || source.remaining <= 0.0) continue;   // not active

    double source_position[3] = { source.x, source.y, source.z };
    for (int iparticle = 0; iparticle < num_particles; ++iparticle)
    {
      double particle_position[3];
      for (int idir = 0; idir < (int)D; ++idir) particle_position[idir] = device_position(iparticle, idir);
      const double w = SPHkernel<D>::kernel(particle_position, source_position,
                                            smearing_radius, smearing_radius_eta);
      if (w <= 0.0) continue;

      // particle a receives (dE W_a / tau) * (tau m_a / sigma_lab_a) = dE W_a m_a / sigma_lab_a,
      // so coverage = sum_a W_a m_a / sigma_lab_a (no tau; same in Cartesian)
      const double vol = device_sph_mass(iparticle, ccake::densities_info::s)
                         / device_hydro_scalar(iparticle, ccake::hydro_info::sigma_lab);
      const double weight = w * vol;
      const int freeze = device_freeze(iparticle);
      source.coverage += weight;
      if (freeze >= 4) source.frozen_weight += weight;
      if (freeze == 4) source.permanently_frozen_weight += weight;
      if (freeze == 5)
      {
        source.frozen_at_birth_particles.push_back(iparticle);
        source.frozen_at_birth_weights.push_back(weight);
      }
    }
  }
}

//==============================================================================
// Deposits onto frozen matter never reach Cooper-Frye. Freeze = 4: lost at once.
// Freeze = 5 (frozen at birth): lost only if it never reheats, so it is kept per
// particle and resolved in the ledger. Particle a gets weight_a / coverage of the
// delivered chunk, as in the deposit.
template <unsigned int D>
void Source<D>::record_frozen_deposit(ICSource& source)
{
  const double delivered = source.dE * delivered_fraction(source.coverage);
  for (std::size_t k = 0; k < source.frozen_at_birth_particles.size(); ++k)
  {
    const double weight = source.frozen_at_birth_weights[k];
    PendingFrozenDeposit& pending = source.pending_frozen_deposits[source.frozen_at_birth_particles[k]];
    pending.energy += delivered * (weight / source.coverage);
    if (source.dump_now) pending.charge_share += weight / source.coverage;
  }
  source.energy_onto_frozen             += delivered * source.frozen_weight / source.coverage;
  source.energy_onto_permanently_frozen += delivered * source.permanently_frozen_weight / source.coverage;
  if (source.dump_now)
  {
    source.frozen_share_at_absorb               = source.frozen_weight / source.coverage;
    source.charge_share_onto_permanently_frozen = source.permanently_frozen_weight / source.coverage;
  }
}

//==============================================================================
// Once per timestep (not per RK substep): give each active source its energy chunk,
// deplete its budget, release its charge when absorbed, and move it in propagate mode.
template <unsigned int D>
void Source<D>::advance_sources()
{
  if (!has_ic_sources || ic_sources.empty()) return;
  ++src_step;   // timestep counter for source-lifetime diagnostics

  const bool renormalize = settingsPtr->source_kernel_renormalize;
  if (renormalize) compute_coverage();

  const double t                 = systemPtr->t;   // current time (tau in hyperbolic, t in Cartesian)
  const double dt                = settingsPtr->dt;
  const double loss_coefficient  = settingsPtr->source_loss_coefficient; // constant loss rate [fm^-2]
  const int    num_deposit_steps = std::max(1, settingsPtr->source_deposit_steps);
  const bool   is_hyperbolic     = (settingsPtr->coordinate_system == "hyperbolic");

  for (ICSource& source : ic_sources)
  {
    source.dE = 0.0;
    source.dump_now = false;
    source.vx = 0.0;
    source.vy = 0.0;
    source.veta = 0.0;

    // not formed yet, or already absorbed
    if (t < source.formation_time || source.remaining <= 0.0) continue;

    if (source.activate_step < 0) source.activate_step = src_step;

    // the fluid has not reached it yet: keep the budget and retry next step
    if (renormalize && source.coverage < settingsPtr->source_min_coverage)
    {
      source.defer_steps++;
      continue;
    }

    // free-streaming velocity: momentum chunk dE*v, and the step in propagate mode
    double vx, vy, veta;
    const double energy_cartesian = std::sqrt(source.px*source.px + source.py*source.py
                                            + source.pz*source.pz + source.mass*source.mass);
    if (is_hyperbolic)
    {
      // Milne velocity p^d/p^tau at the current (tau, eta_s), from the Cartesian momentum
      double eta_s = 0.0;
      if (D == 3) eta_s = source.z;
      const double cosh_eta = std::cosh(eta_s), sinh_eta = std::sinh(eta_s);
      double p_tau = energy_cartesian*cosh_eta - source.pz*sinh_eta; // contravariant Milne energy
      if (p_tau <= 0.0) p_tau = 1e-12;
      const double p_eta = (source.pz*cosh_eta - energy_cartesian*sinh_eta) / t; // uses CURRENT tau
      vx = source.px/p_tau;
      vy = source.py/p_tau;
      veta = p_eta/p_tau;
    }
    else
    {
      double energy = energy_cartesian;
      if (energy <= 0.0) energy = 1e-12;
      vx = source.px/energy;
      vy = source.py/energy;
      veta = source.pz/energy;
    }

    // propagate: constant loss rate; instant: E0/deposit_steps per step at the formation point
    double dE = source.E0 / num_deposit_steps;
    if (settingsPtr->source_propagate) dE = loss_coefficient * dt;
    if (dE > source.remaining) dE = source.remaining;   // last step clamps to 0

    if (renormalize)
    {
      // only the delivered part of the chunk leaves the budget (capped 1/coverage);
      // a negligible tail counts as absorbed so the charge is still released
      source.remaining -= dE * delivered_fraction(source.coverage);
      if (source.remaining < 1.0e-6 * std::max(source.E0, 1.0e-12))
        source.remaining = 0.0;
    }
    else
    {
      source.remaining -= dE;
    }

    source.dE = dE;
    source.vx = vx;
    source.vy = vy;
    source.veta = veta;

    // absorbed: release the charge this step
    if (source.remaining <= 0.0 && !source.charge_dumped)
    {
      source.dump_now      = true;
      source.charge_dumped = true;
      source.absorb_step   = src_step;   // record absorption for lifetime stats
      if (source.coverage > 0.0) source.frozen_share_at_absorb = 0.0;  // record_frozen_deposit sets it
    }

    // before the move: the coverage was computed at this position
    if (renormalize && source.coverage > 0.0) record_frozen_deposit(source);

    if (settingsPtr->source_propagate)
    {
      source.x += vx*dt;
      source.y += vy*dt;
      if (D == 3) source.z += veta*dt;
    }
  }

  // Refresh the ledger (cheap; survives a mid-run crash).
  write_source_ledger();
}

//==============================================================================
// One row per source: identity, kinematics, charges, delivered energy, activation
// and absorption steps, frozen-matter losses. The corona list (energy left or
// absorb_step < 0) and the lifetime totals are derived from it in post-processing.
// Rewritten every step, so a crashed run keeps it up to that step.
template <unsigned int D>
void Source<D>::write_source_ledger()
{
  if (ic_sources.empty()) return;

  std::ofstream ledger(settingsPtr->results_directory.string() + "/source_ledger.dat");
  if (!ledger.is_open()) return;

  ledger << "# per-source ledger at step " << src_step
         << ", tau = " << systemPtr->t << " fm\n";
  ledger << "# source_model " << settingsPtr->source_model << "\n";
  ledger << "# columns: idx type parent PID tau_form x y eta_s E0_GeV "
            "E_delivered_GeV E_left_GeV px_GeV py_GeV pz_GeV mass_GeV "
            "B S Q activate_step absorb_step defer_steps status "
            "E_onto_frozen_GeV frozen_share_at_absorb E_lost_GeV charge_share_lost\n";
  ledger << "# E_onto_frozen / frozen_share_at_absorb: everything that landed on Freeze >= 4 or 5 matter when deposited.\n"
            "# E_lost / charge_share_lost: the part that never reaches Cooper-Frye -- on Freeze = 4 matter, or on\n"
            "# frozen-at-birth particles that are still frozen now (never reheated). Equal to the former without reheating.\n";
  ledger << "# status: 0 dormant(never formed)  1 activated(budget left)  "
            "2 absorbed\n";
  ledger << "# positions are (x, y, eta_s) with the Milne formation time in "
            "hyperbolic mode; momenta are the conserved Cartesian components\n";
  ledger << std::scientific << std::setprecision(8);

  CREATE_VIEW(device_, systemPtr->cabana_particles);
  for (std::size_t i = 0; i < ic_sources.size(); ++i)
  {
    const ICSource& source = ic_sources[i];

    // lost = on Freeze = 4 matter + on frozen-at-birth particles still frozen now
    double energy_lost       = source.energy_onto_permanently_frozen;
    double charge_share_lost = source.charge_share_onto_permanently_frozen;
    for (const auto& entry : source.pending_frozen_deposits)
    {
      const int iparticle = entry.first;
      if (device_freeze(iparticle) != 5) continue;
      energy_lost       += entry.second.energy;
      charge_share_lost += entry.second.charge_share;
    }

    int status = 2;                            // absorbed
    if (source.absorb_step < 0)   status = 1;  // activated, budget left
    if (source.activate_step < 0) status = 0;  // never formed

    const double energy_delivered = source.E0 - source.remaining;
    ledger << i << " " << source.type << " " << source.parent << " "
           << source.PID << " " << source.formation_time << " "
           << source.x << " " << source.y << " " << source.z << " "
           << source.E0*hbarc_GeVfm << " "
           << energy_delivered*hbarc_GeVfm << " "
           << source.remaining*hbarc_GeVfm << " "
           << source.px*hbarc_GeVfm << " "
           << source.py*hbarc_GeVfm << " "
           << source.pz*hbarc_GeVfm << " "
           << source.mass*hbarc_GeVfm << " "
           << source.baryon << " " << source.strange << " " << source.charge << " "
           << source.activate_step << " " << source.absorb_step << " "
           << source.defer_steps << " " << status << " "
           << source.energy_onto_frozen*hbarc_GeVfm << " "
           << source.frozen_share_at_absorb << " "
           << energy_lost*hbarc_GeVfm << " "
           << charge_share_lost << "\n";
  }
  ledger.close();
}

//==============================================================================
// Every RK substep presents the same per-step chunk, so a step deposits exactly dE.
template <unsigned int D>
void Source<D>::build_step_sources()
{
  n_sources = ic_sources.size();
  SRC_VIEW(s_, sources);
  // Shepard renormalization divides by the coverage, capped at max_renormalization_factor
  const bool   renormalize    = settingsPtr->source_kernel_renormalize;
  const double coverage_floor = 1.0 / std::max(1.0, settingsPtr->source_max_renormalization_factor);

  for (std::size_t i = 0; i < n_sources; ++i)
  {
    const ICSource& source = ic_sources[i];
    s_position(i,0) = source.x;
    s_position(i,1) = source.y;
    s_position(i,2) = source.z;
    s_momentum(i,0) = source.dE * source.vx;
    s_momentum(i,1) = source.dE * source.vy;
    s_momentum(i,2) = source.dE * source.veta;
    s_energy(i) = source.dE;               // 0 -> deposit_sources skips it

    // the charge goes in only on the absorption step
    s_baryon_charge(i)      = 0.0;
    s_strangeness_charge(i) = 0.0;
    s_electric_charge(i)    = 0.0;
    if (source.dump_now)
    {
      s_baryon_charge(i)      = source.baryon;
      s_strangeness_charge(i) = source.strange;
      s_electric_charge(i)    = source.charge;
    }

    s_renormalization(i) = 1.0;
    if (renormalize) s_renormalization(i) = std::max(source.coverage, coverage_floor);
  }
}

}// namespace ccake
