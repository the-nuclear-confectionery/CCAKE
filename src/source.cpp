#include "source.h"

namespace ccake
{
// Template instantiations. Source<D> is held by SPHWorkstation<D>, which is
// instantiated for D=1,2,3 (see src/sph_workstation.cpp).
template class Source<1>;
template class Source<2>;
template class Source<3>;

//==============================================================================
template <unsigned int D>
void Source<D>::initialize( std::shared_ptr<Settings> settingsPtr_in,
                            std::shared_ptr<SystemState<D>> systemPtr_in )
{
  settingsPtr = settingsPtr_in;
  systemPtr   = systemPtr_in;
  formatted_output::report("Initializing source terms");

  // check for toymodel/bbmg
  if (settingsPtr->source_model != "toy_model" && settingsPtr->source_model != "BBMG")
  {
    formatted_output::detail("Reading source file");
    read_source_file();
    // IC sources (SMASH particles or AMPT partons) are read once into the
    // master state and deposited through the real-time framework
    // (advance_sources + build_step_sources) for BOTH modes. This presents the
    // same per-step chunk on every RK substep, so the deposit is RK-consistent
    // (the legacy window path fired on different substeps depending on t_src
    // alignment, injecting an inconsistent kick between predictor/corrector that
    // could seed a viscous blow-up).
    //   propagate=true  -> constant-rate energy loss along the trajectory
    //   propagate=false -> instant: budget dumped over deposit_steps at formation
    if (settingsPtr->source_type == "ic" &&
        (settingsPtr->source_model == "SMASH" || settingsPtr->source_model == "AMPT"))
    {
      has_propagating_sources = true;
      formatted_output::detail(settingsPtr->source_model + " sources: "
          + (settingsPtr->source_propagate
             ? "propagate & deplete (constant-rate loss)"
             : "instant deposit (RK-consistent real-time path)"));
    }
  }
  else
  {
    has_jet_sources = true;
    formatted_output::detail("Initializing real-time jet source term");
  }
}

//==============================================================================
template <unsigned int D>
void Source<D>::read_source_file()
{
  if (settingsPtr->source_type == "ic")
  {
    if (settingsPtr->source_model == "AMPT")  read_ampt_file();
    if (settingsPtr->source_model == "SMASH") read_smash_file();
  }
  if (settingsPtr->source_type == "jet")
  {
    if (settingsPtr->source_model == "JEWEL") read_jewel_file();
  }
}

//==============================================================================
template <unsigned int D>
void Source<D>::build_jet_sources(double t, double t_prev)
{
  if (settingsPtr->source_model == "toy_model")
    build_toy_sources(t, t_prev);
  else if (settingsPtr->source_model == "BBMG")
    build_bbmg_sources(t, t_prev);
}

//==============================================================================
template <unsigned int D>
void Source<D>::deposit_sources()
{
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  double t      = systemPtr->t;
  double dt     = settingsPtr->dt;
  double t_prev = t - dt;

  // check for real-time jet source
  if (has_jet_sources)
    build_jet_sources(t, t_prev);
  // Rebuild the per-step source entries for propagating IC sources.
  if (has_propagating_sources)
    build_step_sources(t, t_prev);

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

  // Deposit ALL active, in-window sources in a SINGLE kernel: each SPH particle
  // loops over the sources internally. This is O(1) kernel launches instead of
  // one launch per source -- essential for AMPT (thousands of partons); for
  // the handful of SMASH/BBMG sources it is equivalent. The source AoSoA slices
  // are device-accessible, so the per-source data is read on-device.
  const double normalization        = settingsPtr->source_normalization;
  const bool   is_hyperbolic         = (settingsPtr->coordinate_system == "hyperbolic");
  const bool   deposit_baryon        = settingsPtr->baryon_source;
  const bool   deposit_strangeness   = settingsPtr->strangeness_source;
  const bool   deposit_electric      = settingsPtr->electric_source;
  const int    num_sources           = (int)n_sources;

  // host-side count of active-in-window sources (for the log)
  int num_sources_deposited = 0;
  for (int isrc = 0; isrc < num_sources; ++isrc)
    if (s_energy(isrc) > 0.0 && s_time(isrc) > t_prev && s_time(isrc) <= t) num_sources_deposited++;

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
      if (s_energy(isrc) <= 0.0) continue;                 // exhausted / inactive
      const double source_time = s_time(isrc);
      if (!(source_time > t_prev && source_time <= t)) continue; // not in this window

      double source_position[D];
      for (int idir = 0; idir < D; ++idir) source_position[idir] = s_position(isrc, idir);
      // Split (factorized) deposit kernel: transverse smearing_radius x
      // longitudinal smearing_radius_eta. For D<3 the eta length is ignored
      // (the overload reduces to the isotropic transverse kernel).
      double kernel_value = SPHkernel<D>::kernel(particle_position, source_position,
                                                 smearing_radius, smearing_radius_eta);
      if (kernel_value <= 0.0) continue;
      if (is_hyperbolic) kernel_value /= t;                // Milne volume Jacobian

      const double weight = normalization * kernel_value / dt; // -> Energy/time, Momentum/time
      energy_current += weight * s_energy(isrc);
      for (int idir = 0; idir < D; ++idir) momentum_current[idir] += weight * s_momentum(isrc, idir);
      if (deposit_baryon)      rho_baryon      += s_baryon_charge(isrc)      * kernel_value / dt;
      if (deposit_strangeness) rho_strangeness += s_strangeness_charge(isrc) * kernel_value / dt;
      if (deposit_electric)    rho_electric    += s_electric_charge(isrc)    * kernel_value / dt;
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
void Source<D>::allocate_master_state(std::size_t n)
{
  src_x.resize(n); src_y.resize(n); src_z.resize(n);
  src_px.resize(n); src_py.resize(n); src_pz.resize(n);
  src_mass.resize(n); src_E0.resize(n);
  src_baryon.resize(n); src_strange.resize(n); src_charge.resize(n);
  src_formation_time.resize(n); src_remaining.resize(n);
  src_charge_dumped.assign(n, 0); src_PID.resize(n);
  src_dE.assign(n, 0.0);
  src_vx.assign(n, 0.0); src_vy.assign(n, 0.0); src_veta.assign(n, 0.0);
  src_dump_now.assign(n, 0);
  src_activate_step.assign(n, -1);
  src_absorb_step.assign(n, -1);
}

//==============================================================================
template <unsigned int D>
void Source<D>::register_source(std::size_t i,
    double px, double py, double pz, double E, double mass,
    double x, double y, double z, double t,
    int pid, double baryon_charge, double strangeness_charge, double electric_charge,
    bool seed_master)
{
  const bool is_hyperbolic = (settingsPtr->coordinate_system == "hyperbolic");

  // Keep the conserved CARTESIAN momentum: the Milne velocity is recomputed at
  // the current (tau, eta_s) each step (see advance_sources), not frozen here.
  const double px_cartesian = px, py_cartesian = py, pz_cartesian = pz;
  double E_milne = E;
  transform_to_hyperbolic(is_hyperbolic, px, py, pz, x, y, z, t, E_milne);

  SRC_VIEW(s_, sources);
  s_position(i, 0) = x;
  s_position(i, 1) = y;
  s_position(i, 2) = z;

  s_momentum(i, 0) = px/hbarc_GeVfm;
  s_momentum(i, 1) = py/hbarc_GeVfm;
  s_momentum(i, 2) = pz/hbarc_GeVfm;

  s_energy(i)       = E_milne/hbarc_GeVfm;
  s_total_energy(i) = E_milne/hbarc_GeVfm;
  s_mass(i)         = mass/hbarc_GeVfm;
  s_baryon_charge(i)      = baryon_charge;
  s_strangeness_charge(i) = strangeness_charge;
  s_electric_charge(i)    = electric_charge;
  s_time(i) = t;
  s_PID(i)  = pid;

  if (seed_master)
  {
    // Position in Milne (x, y, eta_s); momentum/mass in CARTESIAN components
    // (fm^-1) so the Milne velocity can be rebuilt each step. The energy budget
    // is the initial Milne energy p^tau (= transformed E).
    src_x[i] = x; src_y[i] = y; src_z[i] = z;          // z = eta_s
    src_px[i] = px_cartesian/hbarc_GeVfm; src_py[i] = py_cartesian/hbarc_GeVfm; src_pz[i] = pz_cartesian/hbarc_GeVfm;
    src_mass[i] = mass/hbarc_GeVfm; src_E0[i] = E_milne/hbarc_GeVfm;
    src_baryon[i] = baryon_charge; src_strange[i] = strangeness_charge; src_charge[i] = electric_charge;
    src_formation_time[i] = t; src_remaining[i] = E_milne/hbarc_GeVfm; src_PID[i] = pid;
  }
}

//==============================================================================
template <unsigned int D>
void Source<D>::log_source_entry(std::size_t i)
{
  SRC_VIEW(s_, sources);
  formatted_output::detail("Source particle " + std::to_string(i) + ": "
    + "Position: (" + std::to_string(s_position(i,0)) + ", " + std::to_string(s_position(i,1)) + ", " + std::to_string(s_position(i,2)) + "), "
    + "Momentum: (" + std::to_string(s_momentum(i,0)) + ", " + std::to_string(s_momentum(i,1)) + ", " + std::to_string(s_momentum(i,2)) + "), "
    + "Energy: " + std::to_string(s_energy(i)) + ", "
    + "Mass: " + std::to_string(s_mass(i)) + ", "
    + "Baryon: " + std::to_string(s_baryon_charge(i)) + ", "
    + "Strangeness: " + std::to_string(s_strangeness_charge(i)) + ", "
    + "Electric: " + std::to_string(s_electric_charge(i)) + ", "
    + "Time: " + std::to_string(s_time(i)) + ", "
    + "PID: " + std::to_string(s_PID(i)));
}

//==============================================================================
template <unsigned int D>
void Source<D>::read_ampt_file()
{
  std::ifstream infile(settingsPtr->source_input_file);
  if (!infile.is_open())
    throw std::runtime_error("Failed to open source file: " + settingsPtr->source_input_file);

  // Skip the first line
  std::string skip_line;
  std::getline(infile, skip_line);

  double px, py, pz, mass, x, y, z, t;
  int pid;
  double dummy1, dummy2, dummy3;

  std::vector<double> px_vec, py_vec, pz_vec, mass_vec, x_vec, y_vec, z_vec, t_vec;
  std::vector<int> pid_vec;

  while (infile >> pid >> px >> py >> pz >> mass >> x >> y >> z >> t >> dummy1 >> dummy2 >> dummy3)
  {
    if (!passes_time_cut(t, z)) continue;
    px_vec.push_back(px);   py_vec.push_back(py);   pz_vec.push_back(pz);
    mass_vec.push_back(mass);
    x_vec.push_back(x);     y_vec.push_back(y);     z_vec.push_back(z);
    t_vec.push_back(t);     pid_vec.push_back(pid);
  }
  infile.close();

  n_sources = px_vec.size();
  sources   = source_array("sources", n_sources);
  allocate_master_state(n_sources);

  for (std::size_t i = 0; i < n_sources; ++i)
  {
    const double E = std::sqrt(px_vec[i]*px_vec[i] + py_vec[i]*py_vec[i]
                             + pz_vec[i]*pz_vec[i] + mass_vec[i]*mass_vec[i]);

    // Charge from the parton flavor code. AMPT's internal flavor convention
    // (as assigned by AMPTGenesis src/AMPT_smearer.cpp) is 1=up, 2=down,
    // 3=strange, 4=charm, 5=bottom, 6=top -- NOT the PDG 1=d, 2=u numbering.
    // baryon_charge = +1/3 for quarks, -1/3 for antiquarks.
    double electric_charge, baryon_charge, strangeness_charge;
    switch (pid_vec[i])
    {
      case 1:   // up quark
        electric_charge = 2.0/3.0;  baryon_charge = 1.0/3.0;  strangeness_charge = 0.0; break;
      case -1:  // anti-up quark
        electric_charge = -2.0/3.0; baryon_charge = -1.0/3.0; strangeness_charge = 0.0; break;
      case 2:   // down quark
        electric_charge = -1.0/3.0; baryon_charge = 1.0/3.0;  strangeness_charge = 0.0; break;
      case -2:  // anti-down quark
        electric_charge = 1.0/3.0;  baryon_charge = -1.0/3.0; strangeness_charge = 0.0; break;
      case 3:   // strange quark
        electric_charge = -1.0/3.0; baryon_charge = 1.0/3.0;  strangeness_charge = -1.0; break;
      case -3:  // anti-strange quark
        electric_charge = 1.0/3.0;  baryon_charge = -1.0/3.0; strangeness_charge = 1.0; break;
      case 4:   // charm quark
        electric_charge = 2.0/3.0;  baryon_charge = 1.0/3.0;  strangeness_charge = 0.0; break;
      case -4:  // anti-charm quark
        electric_charge = -2.0/3.0; baryon_charge = -1.0/3.0; strangeness_charge = 0.0; break;
      case 5:   // bottom quark
        electric_charge = -1.0/3.0; baryon_charge = 1.0/3.0;  strangeness_charge = 0.0; break;
      case -5:  // anti-bottom quark
        electric_charge = 1.0/3.0;  baryon_charge = -1.0/3.0; strangeness_charge = 0.0; break;
      case 6:   // top quark
        electric_charge = 2.0/3.0;  baryon_charge = 1.0/3.0;  strangeness_charge = 0.0; break;
      case -6:  // anti-top quark
        electric_charge = -2.0/3.0; baryon_charge = -1.0/3.0; strangeness_charge = 0.0; break;
      case 21:  // gluon
        electric_charge = 0.0; baryon_charge = 0.0; strangeness_charge = 0.0; break;
      default:
        throw std::runtime_error("Unknown PID in source file");
    }

    register_source(i, px_vec[i], py_vec[i], pz_vec[i], E, mass_vec[i],
                    x_vec[i], y_vec[i], z_vec[i], t_vec[i], pid_vec[i],
                    baryon_charge, strangeness_charge, electric_charge, true);

    //if (i < 10) log_source_entry(i); // AMPT has thousands of partons
    log_source_entry(i);
  }
  formatted_output::detail("Loaded " + std::to_string(n_sources) + " AMPT parton sources");
}

//==============================================================================
template <unsigned int D>
void Source<D>::read_smash_file()
{
  std::ifstream infile(settingsPtr->source_input_file);
  if (!infile.is_open())
    throw std::runtime_error("Failed to open source file: " + settingsPtr->source_input_file);

  std::string line;
  double px, py, pz, mass, x, y, z, t;
  int pid;
  double electric_charge, baryon_charge, strangeness_charge, form_time;
  // Dummy placeholders
  double p0, xsecfac, time_last_coll;
  int idd;
  int ncoll, proc_id_origin, proc_type_origin, pdg_mother1, pdg_mother2;

  std::vector<double> px_vec, py_vec, pz_vec, mass_vec, x_vec, y_vec, z_vec, t_vec;
  std::vector<double> electric_charge_vec, baryon_charge_vec, strangeness_charge_vec;
  std::vector<int> pid_vec;

  while (std::getline(infile, line))
  {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream line_stream(line);
    // Full OSCAR2013Extended particle_lists row. Read electric_charge,
    // baryon_number and strangeness from file so all charges are physical;
    // whether each is actually deposited is gated later by the per-charge
    // enable flags (settingsPtr->{baryon,strangeness,electric}_source).
    if (!(line_stream >> t >> x >> y >> z >> mass >> p0 >> px >> py >> pz >> pid >> idd
              >> electric_charge >> ncoll >> form_time >> xsecfac
              >> proc_id_origin >> proc_type_origin >> time_last_coll
              >> pdg_mother1 >> pdg_mother2 >> baryon_charge >> strangeness_charge))
      continue;

    if (!passes_time_cut(t, z)) continue;
    px_vec.push_back(px);   py_vec.push_back(py);   pz_vec.push_back(pz);
    mass_vec.push_back(mass);
    x_vec.push_back(x);     y_vec.push_back(y);     z_vec.push_back(z);
    t_vec.push_back(t);     pid_vec.push_back(pid);
    electric_charge_vec.push_back(electric_charge);
    baryon_charge_vec.push_back(baryon_charge);
    strangeness_charge_vec.push_back(strangeness_charge);
  }
  infile.close();

  n_sources = px_vec.size();
  sources   = source_array("sources", n_sources);
  allocate_master_state(n_sources);

  for (std::size_t i = 0; i < n_sources; ++i)
  {
    const double E = std::sqrt(px_vec[i]*px_vec[i] + py_vec[i]*py_vec[i]
                             + pz_vec[i]*pz_vec[i] + mass_vec[i]*mass_vec[i]);
    register_source(i, px_vec[i], py_vec[i], pz_vec[i], E, mass_vec[i],
                    x_vec[i], y_vec[i], z_vec[i], t_vec[i], pid_vec[i],
                    baryon_charge_vec[i], strangeness_charge_vec[i], electric_charge_vec[i], true);
    log_source_entry(i);
  }
  formatted_output::detail("Loaded " + std::to_string(n_sources) + " SMASH particle sources");
}

//==============================================================================
template <unsigned int D>
void Source<D>::read_jewel_file()
{
  std::ifstream infile(settingsPtr->source_input_file);
  if (!infile.is_open())
    throw std::runtime_error("Failed to open source file: " + settingsPtr->source_input_file);

  std::string line;

  // JEWEL M line:  M  t  x  y  z  E  px  py  pz   (pos in fm, mom/energy in GeV)
  double t, x, y, z, E, px, py, pz;

  std::vector<double> t_vec, x_vec, y_vec, z_vec;
  std::vector<double> E_vec, px_vec, py_vec, pz_vec;
  std::vector<double> sort_time; // t (cartesian) or tau (hyperbolic), for ordering

  while (std::getline(infile, line))
  {
    if (line.empty() || line[0] == '#') continue;

    std::istringstream line_stream(line);
    std::string tag;
    line_stream >> tag;
    if (!line_stream) continue;
    if (tag != "M") continue;

    if (!(line_stream >> t >> x >> y >> z >> E >> px >> py >> pz)) continue;

    if (!passes_time_cut(t, z)) continue;

    t_vec.push_back(t);   x_vec.push_back(x);   y_vec.push_back(y);   z_vec.push_back(z);
    E_vec.push_back(E);   px_vec.push_back(px); py_vec.push_back(py); pz_vec.push_back(pz);

    // sort key: tau in hyperbolic mode, else t (matches passes_time_cut's start_time)
    const double start_time = (settingsPtr->coordinate_system == "hyperbolic")
                              ? std::sqrt(t*t - z*z) : t;
    sort_time.push_back(start_time);
  }
  infile.close();

  // Order by time once read
  std::vector<std::size_t> sorted_indices(sort_time.size());
  for (std::size_t i = 0; i < sorted_indices.size(); ++i) sorted_indices[i] = i;
  std::sort(sorted_indices.begin(), sorted_indices.end(),
            [&](std::size_t a, std::size_t b){ return sort_time[a] < sort_time[b]; });

  n_sources = sorted_indices.size();
  sources   = source_array("sources", n_sources);

  for (std::size_t i = 0; i < n_sources; ++i)
  {
    const std::size_t source_index = sorted_indices[i];
    // JEWEL M-lines are massless scattering centres: energy comes from the file
    // (not sqrt(p^2+m^2)), charges are discarded, PID is a sentinel, and no
    // propagating master state is seeded.
    register_source(i, px_vec[source_index], py_vec[source_index], pz_vec[source_index],
                    E_vec[source_index], /*mass*/0.0,
                    x_vec[source_index], y_vec[source_index], z_vec[source_index],
                    t_vec[source_index], /*pid*/0,
                    0.0, 0.0, 0.0, /*seed_master*/false);
  }

  formatted_output::detail("Loaded " + std::to_string(n_sources) + " JEWEL M-line sources (sorted by time)");
  for (std::size_t i = 0; i < n_sources; ++i) log_source_entry(i);
}

//==============================================================================
template <unsigned int D>
void Source<D>::build_bbmg_sources(double t, double t_prev)
{
  // Deposit each jet's per-step energy loss (and momentum recoil along its
  // direction) into the fluid. The jet state is already in simulation
  // coordinates, so no transform_to_hyperbolic is applied (unlike the toy
  // model). deposit_sources() then SPH-deposits these via the factorized kernel,
  // and its /dt turns the per-step dE into the energy rate j0_ext expects.
  if (bbmgPtr == nullptr)
  {
    formatted_output::detail("BBMG source requested but no BBMG object was set.");
    return;
  }

  // Refresh the host-side mirror so we can read updated jet positions/e_loss.
  bbmgPtr->copy_device_to_host_BBMG();
  auto& jets = bbmgPtr->jetInfo_host;

  n_sources = jets.size();
  if (n_sources == 0) return;
  sources = source_array("sources", n_sources);
  SRC_VIEW(s_, sources);

  double tau = systemPtr->t;
  for (std::size_t i = 0; i < jets.size(); ++i)
  {
    double dE = jets[i].e_loss;

    // Propagation direction (same construction as BBMG::propagate).
    // direction is a UNIT velocity in the Milne metric:
    //   |direction|^2_Milne = sech^2(delta_y) + tanh^2(delta_y) = 1,
    // so the deposited 4-current J^mu = dE*(1, direction) is null -- energy and
    // momentum match exactly per step (massless transfer, E = |p|_Milne).
    // (Over a curving path direction rotates, so the time-integrated deposit is
    //  timelike; that is the correct physics, not an inconsistency.)
    double delta_y       = (D == 3) ? jets[i].rapidity - jets[i].eta : 0.0;
    double sech_delta_y  = 1.0 / std::cosh(delta_y);
    double direction[3] = { std::cos(jets[i].phi) * sech_delta_y,
                            std::sin(jets[i].phi) * sech_delta_y,
                            (D == 3) ? std::tanh(delta_y) / tau : 0.0 };

    for (int idir = 0; idir < 3; ++idir)
    {
      s_position(i, idir) = (idir < (int)D) ? jets[i].r[idir] : 0.0;
      s_momentum(i, idir) = dE * direction[idir]; // massless: |p|_Milne = E along direction
    }
    s_total_energy(i) = jets[i].total_energy;
    s_energy(i) = dE;
    s_mass(i) = 0.0;
    s_baryon_charge(i) = 0.0;
    s_strangeness_charge(i) = 0.0;
    s_electric_charge(i) = 0.0;
    s_time(i) = t; // deposit this step (deposit_sources window: t_prev < t_src <= t)
    s_PID(i) = jets[i].PID;
  }
}

//==============================================================================
// Advance the propagating sources by ONE timestep: move each source, shed an
// energy chunk (clamped to the remaining budget), deplete the budget, and flag
// the charge dump on the absorption step. Must be called ONCE per timestep (from
// SPHWorkstation::advance_timestep, next to bbmg.propagate()) -- NOT from
// deposit_sources, which fires once per RK substep. build_step_sources() only
// READS the results here. This mirrors the BBMG jetInfo->sources split.
template <unsigned int D>
void Source<D>::advance_sources()
{
  if (!has_propagating_sources) return;
  const std::size_t num_sources = src_x.size();
  if (num_sources == 0) return;
  ++src_step;   // timestep counter for source-lifetime diagnostics

  const double t                = systemPtr->t;   // current time (tau in hyperbolic, t in Cartesian)
  const double dt               = settingsPtr->dt;
  const double loss_coefficient = settingsPtr->source_loss_coefficient; // constant loss rate [fm^-2]
  const bool   is_hyperbolic    = (settingsPtr->coordinate_system == "hyperbolic");

  for (std::size_t i = 0; i < num_sources; ++i)
  {
    src_dE[i] = 0.0;
    src_dump_now[i] = 0;
    src_vx[i] = 0.0; src_vy[i] = 0.0; src_veta[i] = 0.0;

    // Not yet active, or already absorbed: stop moving and stop depositing.
    if (t < src_formation_time[i] || src_remaining[i] <= 0.0) continue;

    // First step this source is active: record activation for lifetime stats.
    if (src_activate_step[i] < 0) src_activate_step[i] = src_step;

    // Free-streaming velocity of the source, stored so build_step_sources can
    // deposit the momentum chunk dE * v along the source direction. It is also
    // the propagation velocity in propagate mode.
    double vx, vy, veta;
    const double energy_cartesian = std::sqrt(src_px[i]*src_px[i] + src_py[i]*src_py[i]
                                            + src_pz[i]*src_pz[i] + src_mass[i]*src_mass[i]);
    if (is_hyperbolic)
    {
      // Milne velocity dx^d/dtau = p^d/p^tau, recomputed at the CURRENT
      // (tau, eta_s) from the conserved Cartesian momentum (NOT frozen at the
      // initial transform). p^tau and p^eta rotate as eta_s evolves and the eta
      // term carries 1/tau(current). For D<3 (boost-invariant) eta_s stays 0 and
      // this reduces to pure transverse Cartesian motion.
      const double eta_s = (D == 3) ? src_z[i] : 0.0;
      const double cosh_eta = std::cosh(eta_s), sinh_eta = std::sinh(eta_s);
      double p_tau = energy_cartesian*cosh_eta - src_pz[i]*sinh_eta; // contravariant Milne energy
      if (p_tau <= 0.0) p_tau = 1e-12;
      const double p_eta = (src_pz[i]*cosh_eta - energy_cartesian*sinh_eta) / t; // uses CURRENT tau
      vx = src_px[i]/p_tau; vy = src_py[i]/p_tau; veta = p_eta/p_tau;
    }
    else
    {
      // Cartesian free-streaming: dx/dt = p/E (v/c), no tau/eta geometry.
      double energy = (energy_cartesian > 0.0) ? energy_cartesian : 1e-12;
      vx = src_px[i]/energy; vy = src_py[i]/energy; veta = src_pz[i]/energy;
    }

    // Propagate the source one step ONLY in propagate mode -- streaming the
    // deposit across cells along the path keeps the local shear sigma = grad(u)
    // bounded. In instant mode the source stays put at its formation position:
    // a true instant deposit, spread over deposit_steps in TIME only, not space.
    if (settingsPtr->source_propagate)
    {
      src_x[i] += vx*dt;
      src_y[i] += vy*dt;
      if (D == 3) src_z[i] += veta*dt;
    }

    double dE;
    if (settingsPtr->source_propagate)
    {
      // CONSTANT loss rate dE = loss_coefficient * dt
      dE = loss_coefficient * dt;
      if (dE > src_remaining[i]) dE = src_remaining[i];
    }
    else
    {
      // Instant burst: lay the full budget down as E0/N over N steps at the
      // FIXED formation position (no propagation), so it is spread in time only.
      // Charge releases on the last step (remaining -> 0). deposit_steps=1 dumps
      // the whole budget in a single step.
      const int num_deposit_steps = (settingsPtr->source_deposit_steps > 1)
                                    ? settingsPtr->source_deposit_steps : 1;
      dE = src_E0[i] / num_deposit_steps;
      if (dE > src_remaining[i]) dE = src_remaining[i]; // last step clamps to 0
    }
    src_remaining[i] -= dE;

    src_dE[i] = dE;
    src_vx[i] = vx; src_vy[i] = vy; src_veta[i] = veta;

    // Flag the one-shot charge release on the absorption step.
    if (src_remaining[i] <= 0.0 && !src_charge_dumped[i])
    {
      src_dump_now[i]      = 1;
      src_charge_dumped[i] = 1;
      src_absorb_step[i]   = src_step;   // record absorption for lifetime stats
    }
  }

  // Refresh the running diagnostics file (cheap; survives a mid-run crash).
  write_source_report();
}

//==============================================================================
// Lifetime / escape diagnostics for the propagating sources. Reports the average
// number of timesteps from activation to absorption, and flags sources that
// escape with energy budget (and hence baryon charge) left undeposited. Written
// to <results>/smash_source_report.txt and echoed to stdout.
template <unsigned int D>
void Source<D>::write_source_report()
{
  const std::size_t num_sources = src_remaining.size();
  if (num_sources == 0) return;

  long num_activated = 0, num_absorbed = 0, num_escaped = 0, num_dormant = 0;
  long total_steps = 0, min_steps = -1, max_steps = -1;
  double escaped_energy = 0.0, escaped_baryon = 0.0, escaped_baryon_abs = 0.0, total_baryon_abs = 0.0;
  double dormant_baryon_abs = 0.0;   // |B| of sources that never formed

  for (std::size_t i = 0; i < num_sources; ++i)
  {
    total_baryon_abs += std::abs(src_baryon[i]);
    if (src_activate_step[i] < 0) { ++num_dormant; dormant_baryon_abs += std::abs(src_baryon[i]); continue; } // formation tau never reached
    ++num_activated;
    if (src_absorb_step[i] >= 0)
    {
      ++num_absorbed;
      const int steps = src_absorb_step[i] - src_activate_step[i];
      total_steps += steps;
      if (min_steps < 0 || steps < min_steps) min_steps = steps;
      if (steps > max_steps) max_steps = steps;
    }
    else
    {
      // Active but never absorbed -> escaped with budget (and charge) left.
      ++num_escaped;
      escaped_energy     += src_remaining[i];
      escaped_baryon     += src_baryon[i];
      escaped_baryon_abs += std::abs(src_baryon[i]);
    }
  }
  const double avg_steps_to_absorb = (num_absorbed > 0) ? double(total_steps)/num_absorbed : 0.0;
  const double avg_fm_to_absorb    = avg_steps_to_absorb * settingsPtr->dt;

  std::ostringstream report;
  report << "# SMASH propagating-source lifetime report (step " << src_step
         << ", tau = " << systemPtr->t << " fm)\n";
  report << "total_sources          " << num_sources << "\n";
  report << "dormant(never_formed)  " << num_dormant << "\n";
  report << "activated              " << num_activated << "\n";
  report << "absorbed               " << num_absorbed << "\n";
  report << "escaped(budget_left)   " << num_escaped << "\n";
  report << "avg_steps_to_absorb    " << avg_steps_to_absorb << "   (" << avg_fm_to_absorb << " fm)\n";
  report << "min_steps_to_absorb    " << min_steps << "\n";
  report << "max_steps_to_absorb    " << max_steps << "\n";
  report << "escaped_energy_left    " << escaped_energy << " fm^-1\n";
  report << "escaped_baryon_net     " << escaped_baryon << "\n";
  report << "escaped_baryon_abs     " << escaped_baryon_abs << "  (activated but not absorbed)\n";
  report << "dormant_baryon_abs     " << dormant_baryon_abs << "  (never formed)\n";
  const double baryon_never_deposited = escaped_baryon_abs + dormant_baryon_abs; // escaped + dormant
  report << "baryon_never_deposited " << baryon_never_deposited << " of " << total_baryon_abs
         << "  (" << (total_baryon_abs > 0 ? 100.0*baryon_never_deposited/total_baryon_abs : 0.0)
         << "% of |B| -- only absorbed sources dump their charge)\n";

  const std::string path = settingsPtr->results_directory.string()
                           + "/smash_source_report.txt";
  std::ofstream report_file(path);
  if (report_file.is_open()) { report_file << report.str(); report_file.close(); }
}

//==============================================================================
// Per-substep builder: write the CURRENT master state into the `sources` AoSoA
// that deposit_sources() deposits. Read-only w.r.t. the master state -- all
// advancing/depleting happens once per timestep in advance_sources(). Presenting
// the same per-step chunk (src_dE) across every RK substep makes the RK
// integration deposit exactly src_dE over the timestep (same as BBMG).
template <unsigned int D>
void Source<D>::build_step_sources(double t, double t_prev)
{
  const std::size_t num_sources = src_x.size();
  if (num_sources == 0) return;
  n_sources = num_sources;
  // `sources` was already sized num_sources in the reader; overwrite in place.
  SRC_VIEW(s_, sources);

  for (std::size_t i = 0; i < num_sources; ++i)
  {
    s_position(i,0) = src_x[i];
    s_position(i,1) = src_y[i];
    s_position(i,2) = src_z[i];
    const double dE = src_dE[i];
    s_momentum(i,0) = dE * src_vx[i];
    s_momentum(i,1) = dE * src_vy[i];
    s_momentum(i,2) = dE * src_veta[i];
    s_energy(i) = dE;                      // 0 -> deposit_sources skips it
    s_total_energy(i) = src_remaining[i];
    s_mass(i) = src_mass[i];
    s_time(i) = t;                         // deposit this step (window t_prev < t <= t)
    s_PID(i) = src_PID[i];

    // Charge only on the absorption step (else zero).
    if (src_dump_now[i])
    {
      s_baryon_charge(i)      = src_baryon[i];
      s_strangeness_charge(i) = src_strange[i];
      s_electric_charge(i)    = src_charge[i];
    }
    else
    {
      s_baryon_charge(i) = 0.0; s_strangeness_charge(i) = 0.0; s_electric_charge(i) = 0.0;
    }
  }
}

//==============================================================================
template <unsigned int D>
void Source<D>::build_toy_sources(double t, double t_prev)
{
  // Build sources for the CURRENT timestep window only, so deposit_sources()
  // logic (t_prev < t_src <= t) works unchanged.
  const bool is_hyperbolic = (settingsPtr->coordinate_system == "hyperbolic");

  n_sources = 2;
  sources   = source_array("sources", n_sources);
  SRC_VIEW(s_, sources);

  // Jet starts at x = 0 when t = t0
  double x1 = +(t - settingsPtr->t0), y1 = 0.0, z1 = 0.0;
  double x2 = -(t - settingsPtr->t0), y2 = 0.0, z2 = 0.0;

  // Store time so it is deposited THIS step
  double t_src = t;

  // Toy parameters
  double source_energy = 300.*GeVtofm;

  // Momentum direction: along +/- x
  double px1 = +source_energy, py1 = 0.0, pz1 = 0.0;
  double px2 = -source_energy, py2 = 0.0, pz2 = 0.0;

  double mass = 0.0;
  double E1 = source_energy, E2 = source_energy;

  // Charges ignored for toy model
  double baryon_charge = 0.0, strangeness_charge = 0.0, electric_charge = 0.0;

  // Transform if running in hyperbolic coordinates
  transform_to_hyperbolic(is_hyperbolic, px1, py1, pz1, x1, y1, z1, t_src, E1);

  double t_src2 = t;
  transform_to_hyperbolic(is_hyperbolic, px2, py2, pz2, x2, y2, z2, t_src2, E2);

  // ---- source 0
  s_position(0,0) = x1; s_position(0,1) = y1; s_position(0,2) = z1;
  s_momentum(0,0) = px1; s_momentum(0,1) = py1; s_momentum(0,2) = pz1;
  s_energy(0) = E1;
  s_mass(0) = mass;
  s_baryon_charge(0) = baryon_charge;
  s_strangeness_charge(0) = strangeness_charge;
  s_electric_charge(0) = electric_charge;
  s_time(0) = t_src;
  s_PID(0) = 0;

  // ---- source 1
  s_position(1,0) = x2; s_position(1,1) = y2; s_position(1,2) = z2;
  s_momentum(1,0) = px2; s_momentum(1,1) = py2; s_momentum(1,2) = pz2;
  s_energy(1) = E2;
  s_mass(1) = mass;
  s_baryon_charge(1) = baryon_charge;
  s_strangeness_charge(1) = strangeness_charge;
  s_electric_charge(1) = electric_charge;
  s_time(1) = t_src2;
  s_PID(1) = 0;
}

}// namespace ccake
