#ifndef SOURCE_H
#define SOURCE_H

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

#include <Cabana_Core.hpp>

#include "constants.h"
//#include "eom.h"
#include "eom_default.h"
#include "eom_cartesian.h"
#include "formatted_output.h"
#include "kernel.h"
#include "settings.h"
#include "system_state.h"
#include "constants.h"
#include "thermodynamic_info.h"



namespace ccake
{

#ifndef CONCAT
  #define CONCAT_(a, b) a##b
  #define CONCAT(a, b) CONCAT_(a, b)
#endif
  using SRC = Cabana::MemberTypes<double[3], // position
                                  double[3], // momentum
                                  double,    // energy
                                  double,    // mass
                                  double,    // baryon charge
                                  double,    // strangeness charge
                                  double,    // electric charge
                                  double,    // time
                                  int>;       // PID
                                  
  namespace SRC_enum{
    enum SRC_members
    {
      position = 0,
      momentum,
      energy,
      mass,
      baryon_charge,
      strangeness_charge,
      electric_charge,
      time,
      PID
    };
  }

  using source_array = Cabana::AoSoA<SRC, DeviceType,VECTOR_LENGTH>;

  #define SRC_VIEW(prefix, SRC_aosoa) \
  auto CONCAT(prefix, position) = Cabana::slice<SRC_enum::position>(SRC_aosoa); \
  auto CONCAT(prefix, momentum) = Cabana::slice<SRC_enum::momentum>(SRC_aosoa); \
  auto CONCAT(prefix, energy) = Cabana::slice<SRC_enum::energy>(SRC_aosoa); \
  auto CONCAT(prefix, mass) = Cabana::slice<SRC_enum::mass>(SRC_aosoa); \
  auto CONCAT(prefix, baryon_charge) = Cabana::slice<SRC_enum::baryon_charge>(SRC_aosoa); \
  auto CONCAT(prefix, strangeness_charge) = Cabana::slice<SRC_enum::strangeness_charge>(SRC_aosoa); \
  auto CONCAT(prefix, electric_charge) = Cabana::slice<SRC_enum::electric_charge>(SRC_aosoa); \
  auto CONCAT(prefix, time) = Cabana::slice<SRC_enum::time>(SRC_aosoa); \
  auto CONCAT(prefix, PID) = Cabana::slice<SRC_enum::PID>(SRC_aosoa); \



template <unsigned int D>
class Source
{
private:

  bool real_time_jet_source = false; ///< Flag to indicate if real-time jet source is used

  double total_energy_transferred; ///< Total energy transferred by the source term

  double total_baryon_charge_transferred; ///< Total baryon charge transferred by the source term
  double total_strangeness_charge_transferred; ///< Total strangeness charge transferred by the source term
  double total_electric_charge_transferred; ///< Total electric charge transferred by the source term

  double n_sources; ///< Number of sources in the simulation

  source_array sources; ///< Array of source terms for the particles

  std::shared_ptr<Settings> settingsPtr    = nullptr;     ///< Object containing settings parsed from input file
  std::shared_ptr<SystemState<D>> systemPtr= nullptr; ///< Object containing the SPH System (linked list, particles, etc.)
public:
  /// @brief Constructor for the Source class.
  Source( std::shared_ptr<Settings> settingsPtr_in,
           std::shared_ptr<SystemState<D>> systemPtr_in )
           {
            initialize(settingsPtr_in, systemPtr_in);
           }

  ~Source(){};


  void initialize( std::shared_ptr<Settings> settingsPtr_in,
                              std::shared_ptr<SystemState<D>> systemPtr_in )
  {
    settingsPtr = settingsPtr_in;
    systemPtr = systemPtr_in;
    std::cout << "Initializing source terms" << std::endl;
    //check for toymodel/bbmg 
    if (settingsPtr->source_model != "toy_model" || settingsPtr->source_model != "BBMG")
    {
      //read source file
      std::cout << "Reading source file" << std::endl;
      read_source_file();
    }
    else{
      //initialize source array
      real_time_jet_source = true;
      std::cout << "Initializing real-time jet source term" << std::endl;
    }

  }

  void build_real_time_jet(double t, double t_prev)
  {
    if(settingsPtr->source_model == "toy_model")
    {
      calculate_source_toymodel(t, t_prev);
    }
    else if(settingsPtr->source_model == "BBMG")
    {
      calculate_source_bbmg(t, t_prev);
    }
  }

  void read_source_file()
  {


    if (settingsPtr->source_type == "ic")
    {
      if (settingsPtr->source_model == "AMPT")
      {
        read_source_file_ampt();
     }
      if (settingsPtr->source_model == "SMASH")
      {
        read_source_file_smash();
      }
    }
    if (settingsPtr->source_type == "jet")
    {
     if (settingsPtr->source_model == "JEWEL")
      {
        read_source_file_jewel();
      }
    }
    
  }


  void add_source()
  {
    CREATE_VIEW(device_, systemPtr->cabana_particles);
    auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

    double hT = settingsPtr->hT;
    double t = systemPtr->t;
    double dt = settingsPtr->dt;
    double t_prev = t - dt;
    std::cout << "Adding source terms at time " << t << std::endl;
    std::cout << "Current time: " << t << std::endl;
    std::cout << "Previous time: " << t_prev << std::endl;

    //check for real_time_jet_source
    if (real_time_jet_source)
    {
      build_real_time_jet(t, t_prev);
    }
    
    double smearing_radius = settingsPtr->smearing_radius;
    double n_particles = systemPtr->cabana_particles.size();
    //rest source terms
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
    std::cout << "Adding source terms" << std::endl;
    int n_sources_added = 0;
    for (std::size_t is = 0; is < n_sources; ++is)
    {
      double t_src = s_time(is);
      double dt = t - t_prev; // Time step size, used to determine if source is in the current timestep window
      // Only act if source is within this timestep window
      if (t_src > t_prev && t_src <= t)
      {
        n_sources_added++;
        //double weight = (t - t_src) / dt; // Linear interpolation weight (can reverse if desired)
        double weight = 1.0; // Full weight for simplicity
        double src_pos[D];
        double src_mom[D];
        for (int d = 0; d < D; ++d)
        {
          src_pos[d] = s_position(is, d) * weight;
          src_mom[d] = s_momentum(is, d) * weight;
        }

        double src_energy = s_energy(is) * weight;
        double src_baryon = s_baryon_charge(is) * weight;
        double src_strangeness = s_strangeness_charge(is) * weight;
        double src_charge = s_electric_charge(is) * weight;
        double normalization = settingsPtr->source_normalization; // IC normalization factor
        auto deposit_source = KOKKOS_LAMBDA(const int is, const int ia)
        {
          double r[D];
          for (int idir = 0; idir < D; ++idir)
            r[idir] = device_position.access(is, ia, idir);

          double dist = SPHkernel<D>::distance(r, src_pos);
          double kern = SPHkernel<D>::kernel(dist, smearing_radius);

          if (dist < 2.0 * smearing_radius) // optional cutoff
          {
                      /// \todo Do we need a 1/tau in hyperbolic coordinates here?
            if(settingsPtr->coordinate_system == "hyperbolic")
              kern /= t;
            // Reset accumulation if needed
            device_hydro_scalar.access(is, ia, ccake::hydro_info::j0_ext) += normalization * src_energy * kern /dt; // right units of Energy/time



            for (int idir = 0; idir < D; ++idir)
            {
              device_hydro_vector.access(is, ia, ccake::hydro_info::j_ext, idir) +=normalization * src_mom[idir] * kern / dt; // right units of Momentum/time
            }
            // Optionally deposit charges
            // if (settingsPtr->baryon_source) device_hydro_scalar.access(is, ia, ccake::hydro_info::rho_B_ext) += src_baryon * kern ;
            // if (settingsPtr->strangeness_source) device_hydro_scalar.access(is, ia, ccake::hydro_info::rho_S_ext) += src_strangeness * kern ;
            // if (settingsPtr->electric_source) device_hydro_scalar.access(is, ia, ccake::hydro_info::rho_Q_ext) += src_charge * kern * electric_charge;
          }
        };

        Cabana::simd_parallel_for(simd_policy, deposit_source, "deposit_source");
        Kokkos::fence();
      }
    }

    std::cout << "Added source terms for " << n_sources_added << " source particles out of " << n_sources << std::endl;
  }

  void transform_to_hyperbolic(
      double& px, double& py, double& pz,
      double& x,  double& y,  double& z,
      double& t,
      double& E
  ) const
  {
    if (settingsPtr->coordinate_system != "hyperbolic")
      return;

    if (t*t < z*z)
      return;

    double tau = std::sqrt(t*t - z*z);
    double eta_s = 0.5 * std::log((t + z) / (t - z));


    double p_tau =  E * std::cosh(eta_s) - pz * std::sinh(eta_s);
    double p_eta = (-E * std::sinh(eta_s) + pz * std::cosh(eta_s))/tau;

    // Overwrite momentum & coordinates
    pz = p_eta;
    E  = p_tau;
    t  = tau;
    z  = eta_s;
  }

  void read_source_file_ampt()
  {
    std::ifstream infile(settingsPtr->source_input_file);
    if (!infile.is_open())
    {
      throw std::runtime_error("Failed to open source file: " + settingsPtr->source_input_file);
    }

    // Skip the first line
    std::string skip_line;
    std::getline(infile, skip_line);

    // Now start reading data
    double px, py, pz, mass, x, y, z, t;
    int pid;
    double dummy1, dummy2, dummy3;

    std::vector<double> px_vec, py_vec, pz_vec, mass_vec, x_vec, y_vec, z_vec, t_vec;
    std::vector<int> pid_vec;

    while (infile >> pid >> px >> py >> pz >> mass >> x >> y >> z >> t >> dummy1 >> dummy2 >> dummy3)
    {
      // Check if time is greater than t0
      double x0 = t;
      if (settingsPtr->coordinate_system == "hyperbolic")
      {
        // only add particles in the light cone
        if (t*t < z*z)
        {
          // Skip this entry, as it is outside the light cone in hyperbolic coordinates
          std::cout << "Skipping source entry with t^2 < z^2 at t = " << t << " and z = " << z << std::endl;
          continue;
        }
        double tau = std::sqrt(t*t - z*z);
        double eta_s = 0.5 * std::log((t + z) / (t - z));
        x0 = tau; // Use tau for the comparison in hyperbolic coordinates

      }
      if (x0 > settingsPtr->t0)
      {
        px_vec.push_back(px);
        py_vec.push_back(py);
        pz_vec.push_back(pz);
        mass_vec.push_back(mass);
        x_vec.push_back(x);
        y_vec.push_back(y);
        z_vec.push_back(z);
        t_vec.push_back(t);
        pid_vec.push_back(pid);
      }
    }

    infile.close();

    // Allocate AoSoA
    n_sources = px_vec.size();
    sources = source_array("sources", n_sources);

    // Grab views using your macro
    SRC_VIEW(s_, sources);  // e.g., s_position, s_momentum, etc.

    for (std::size_t i = 0; i < n_sources; ++i)
    {


      double px = px_vec[i];
      double py = py_vec[i];
      double pz = pz_vec[i];
      double mass = mass_vec[i];
      double x = x_vec[i];
      double y = y_vec[i];
      double z = z_vec[i];
      double t = t_vec[i];

      double E = std::sqrt(px*px + py*py + pz*pz + mass*mass);
      transform_to_hyperbolic(px, py, pz, x, y, z, t, E);

      //check the charge against partons PID
      double electric_charge, baryon_charge, strangeness_charge;
      switch (pid_vec[i])
      {
        case 1: // u quark
          electric_charge = 2.0 / 3.0;
          baryon_charge = 1.0 / 3.0;
          strangeness_charge = 0.0;
          break;
        case -1: // d quark
          electric_charge = -1.0 / 3.0;
          baryon_charge = 1.0 / 3.0;
          strangeness_charge = 0.0;
          break;
        case 2: // s quark
          electric_charge = -1.0 / 3.0;
          baryon_charge = 1.0 / 3.0;
          strangeness_charge = -1.0;
          break;
        case -2: // anti-s quark
          electric_charge = 1.0 / 3.0;
          baryon_charge = -1.0 / 3.0;
          strangeness_charge = 1.0;
          break;
        case 3: // c quark
          electric_charge = 2.0 / 3.0;
          baryon_charge = 1.0 / 3.0;
          strangeness_charge = 0.0;
          break;
        case -3: // anti-c quark
          electric_charge = -2.0 / 3.0;
          baryon_charge = -1.0 / 3.0;
          strangeness_charge = 0.0;
          break;
        case 21: // gluon
          electric_charge = 0.0;
          baryon_charge = 0.0;
          strangeness_charge = 0.0;
          break;
        default:
          throw std::runtime_error("Unknown PID in source file");
      }

      s_position(i, 0) = x;
      s_position(i, 1) = y;
      s_position(i, 2) = z;

      s_momentum(i, 0) = px/hbarc_GeVfm; // convert to GeV/fm (if needed, otherwise just use 1.0)
      s_momentum(i, 1) = py/hbarc_GeVfm; // convert to GeV/fm (if needed, otherwise just use 1.0)
      s_momentum(i, 2) = pz/hbarc_GeVfm; // convert to GeV/fm (if needed, otherwise just use 1.0)

      s_energy(i) = E/hbarc_GeVfm;
      s_mass(i) = mass;
      s_baryon_charge(i) = baryon_charge;
      s_strangeness_charge(i) = strangeness_charge;
      s_electric_charge(i) = electric_charge;
      s_time(i) = t;
      s_PID(i) = pid_vec[i];

      //print all info
      std::cout << "Source particle " << i << ": "
                << "Position: (" << x << ", " << y << ", " << z << "), "
                << "Momentum: (" << px << ", " << py << ", " << pz << "), "
                << "Energy: " << E << ", "
                << "Mass: " << mass_vec[i] << ", "
                << "Baryon Charge: " << baryon_charge << ", "
                << "Strangeness Charge: " << strangeness_charge << ", "
                << "Electric Charge: " << electric_charge << ", "
                << "Time: " << t  << ", "
                << "PID: "  << pid_vec[i]  << std::endl;
    }
  }

  void read_source_file_smash()
  {
    std::ifstream infile(settingsPtr->source_input_file);
    if (!infile.is_open())
    {
      throw std::runtime_error("Failed to open source file: " + settingsPtr->source_input_file);
    }

    // Read and process each line, skipping comments
    std::string line;
    double dummy;
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
      std::istringstream iss(line);
      //if (!(iss >> t >> x >> y >> z >> mass >> p0 >> px >> py >> pz >> pid >> idd
      //          >> electric_charge >> ncoll >> form_time >> xsecfac
      //          >> proc_id_origin >> proc_type_origin >> time_last_coll
      //          >> pdg_mother1 >> pdg_mother2 >> baryon_charge >> strangeness_charge))
      if (!(iss >> t >> x >> y >> z >> mass >> p0 >> px >> py >> pz >> pid >> idd
          >> electric_charge ))
      {
        continue;
      }

      double x0 = t;
      if (settingsPtr->coordinate_system == "hyperbolic")
      {
        if (t*t < z*z)
        {
          std::cout << "Skipping source entry with t^2 < z^2 at t = " << t << " and z = " << z << std::endl;
          continue;
        }
        double tau = std::sqrt(t*t - z*z);
        double eta_s = 0.5 * std::log((t + z) / (t - z));
        x0 = tau;
      }
      if (x0 > settingsPtr->t0)
      {
        baryon_charge = 0.0;
        strangeness_charge = 0.0;
        px_vec.push_back(px);
        py_vec.push_back(py);
        pz_vec.push_back(pz);
        mass_vec.push_back(mass);
        x_vec.push_back(x);
        y_vec.push_back(y);
        z_vec.push_back(z);
        t_vec.push_back(t);
        pid_vec.push_back(pid);
        electric_charge_vec.push_back(electric_charge);
        baryon_charge_vec.push_back(baryon_charge);
        strangeness_charge_vec.push_back(strangeness_charge);
      }
    }

    infile.close();

    n_sources = px_vec.size();
    sources = source_array("sources", n_sources);

    SRC_VIEW(s_, sources);

    for (std::size_t i = 0; i < n_sources; ++i)
    {
      double px = px_vec[i];
      double py = py_vec[i];
      double pz = pz_vec[i];
      double mass = mass_vec[i];
      double x = x_vec[i];
      double y = y_vec[i];
      double z = z_vec[i];
      double t = t_vec[i];
      double E = std::sqrt(px*px + py*py + pz*pz + mass*mass);

      transform_to_hyperbolic(px, py, pz, x, y, z, t, E);

      s_position(i, 0) = x;
      s_position(i, 1) = y;
      s_position(i, 2) = z;

      s_momentum(i, 0) = px/hbarc_GeVfm;
      s_momentum(i, 1) = py/hbarc_GeVfm;
      s_momentum(i, 2) = pz/hbarc_GeVfm;

      s_energy(i) = E/hbarc_GeVfm;
      s_mass(i) = mass/hbarc_GeVfm;
      s_baryon_charge(i) = baryon_charge_vec[i];
      s_strangeness_charge(i) = strangeness_charge_vec[i];
      s_electric_charge(i) = electric_charge_vec[i];
      s_time(i) = t;
      s_PID(i) = pid_vec[i];
      //print the loaded sources
      std::cout << "Source particle " << i << ": "
                << "Position: (" << x << ", " << y << ", " << z << "), "
                << "Momentum: (" << px << ", " << py << ", " << pz << "), "
                << "Energy: " << E << ", "
                << "Mass: " << mass << ", "
                << "Baryon Charge: " << baryon_charge << ", "
                << "Strangeness Charge: " << strangeness_charge << ", "
                << "Electric Charge: " << electric_charge << ", "
                << "Time: " << t  << ", "
                << "PID: "  << pid_vec[i]  << std::endl;
    }
  }


  void calculate_source_bbmg(double t, double t_prev)
  {
    // Placeholder for BBMG source calculation
    // Implement BBMG source term logic here
    std::cout << "BBMG source calculation not yet implemented." << std::endl;
  }

  void calculate_source_toymodel(double t, double t_prev)
  {
    // Build sources for the CURRENT timestep window only.
    // This way add_source() logic (t_prev < t_src <= t) works unchanged.

    n_sources = 2.0;
    sources = source_array("sources", n_sources);

    SRC_VIEW(s_, sources);

    // Jet starts at x = 0 when t = t0
    double x1 = +(t - settingsPtr->t0);
    double y1 = 0.0;
    double z1 = 0.0;

    double x2 = -(t - settingsPtr->t0);
    double y2 = 0.0;
    double z2 = 0.0;

    // Store time so it is deposited THIS step
    double t_src = t;

    // Toy parameters
    double source_energy = 300.*GeVtofm;

    // Momentum direction: along ±x
    double px1 = +source_energy;
    double py1 = 0.0;
    double pz1 = 0.0;

    double px2 = -source_energy;
    double py2 = 0.0;
    double pz2 = 0.0;

    double mass = 0.0;
    double E1 = source_energy;
    double E2 = source_energy;

    // Charges ignored for toy model
    double baryon_charge = 0.0;
    double strangeness_charge = 0.0;
    double electric_charge = 0.0;

    // Transform if running in hyperbolic coordinates
    transform_to_hyperbolic(px1, py1, pz1, x1, y1, z1, t_src, E1);

    double t_src2 = t;
    transform_to_hyperbolic(px2, py2, pz2, x2, y2, z2, t_src2, E2);

    // ---- source 0
    s_position(0,0) = x1;
    s_position(0,1) = y1;
    s_position(0,2) = z1;

    s_momentum(0,0) = px1;
    s_momentum(0,1) = py1;
    s_momentum(0,2) = pz1;

    s_energy(0) = E1;
    s_mass(0) = mass;
    s_baryon_charge(0) = baryon_charge;
    s_strangeness_charge(0) = strangeness_charge;
    s_electric_charge(0) = electric_charge;
    s_time(0) = t_src;
    s_PID(0) = 0;

    // ---- source 1
    s_position(1,0) = x2;
    s_position(1,1) = y2;
    s_position(1,2) = z2;

    s_momentum(1,0) = px2;
    s_momentum(1,1) = py2;
    s_momentum(1,2) = pz2;

    s_energy(1) = E2;
    s_mass(1) = mass;
    s_baryon_charge(1) = baryon_charge;
    s_strangeness_charge(1) = strangeness_charge;
    s_electric_charge(1) = electric_charge;
    s_time(1) = t_src2;
    s_PID(1) = 0;
  }


  void read_source_file_jewel()
  {
    std::ifstream infile(settingsPtr->source_input_file);
    if (!infile.is_open())
    {
      throw std::runtime_error("Failed to open source file: " + settingsPtr->source_input_file);
    }

    std::string line;

    // JEWEL M line:  M  t  x  y  z  E  px  py  pz   (pos in fm, mom/energy in GeV)
    double t, x, y, z, E, px, py, pz;

    std::vector<double> t_vec, x_vec, y_vec, z_vec;
    std::vector<double> E_vec, px_vec, py_vec, pz_vec;
    std::vector<double> sort_time; // t (cartesian) or tau (hyperbolic), for ordering

    while (std::getline(infile, line))
    {
      if (line.empty() || line[0] == '#') continue;

      std::istringstream iss(line);
      std::string tag;
      iss >> tag;
      if (!iss) continue;

      if (tag != "M") continue;

      if (!(iss >> t >> x >> y >> z >> E >> px >> py >> pz))
        continue;

      // Decide time for cuts/order (consistent with your other readers)
      double x0 = t;
      if (settingsPtr->coordinate_system == "hyperbolic")
      {
        if (t*t < z*z)
        {
          std::cout << "Skipping source entry with t^2 < z^2 at t = " << t << " and z = " << z << std::endl;
          continue;
        }
        double tau = std::sqrt(t*t - z*z);
        x0 = tau;
      }

      if (x0 > settingsPtr->t0)
      {
        t_vec.push_back(t);
        x_vec.push_back(x);
        y_vec.push_back(y);
        z_vec.push_back(z);

        E_vec.push_back(E);
        px_vec.push_back(px);
        py_vec.push_back(py);
        pz_vec.push_back(pz);

        sort_time.push_back(x0);
      }
    }

    infile.close();

    // Order by time once read
    std::vector<std::size_t> idx(sort_time.size());
    for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;

    std::sort(idx.begin(), idx.end(),
              [&](std::size_t a, std::size_t b){ return sort_time[a] < sort_time[b]; });

    // Allocate AoSoA
    n_sources = idx.size();
    sources = source_array("sources", n_sources);

    SRC_VIEW(s_, sources);

    for (std::size_t i = 0; i < n_sources-4; ++i)
    {
      std::size_t j = idx[i];

      double t = t_vec[j];
      double x = x_vec[j];
      double y = y_vec[j];
      double z = z_vec[j];

      double E  = E_vec[j];
      double px = px_vec[j];
      double py = py_vec[j];
      double pz = pz_vec[j];

      double mass = 0.0;

      // If you want to reuse the same transform helper:
      transform_to_hyperbolic(px, py, pz, x, y, z, t, E);

      s_position(i, 0) = x;
      s_position(i, 1) = y;
      s_position(i, 2) = z;

      s_momentum(i, 0) = px/hbarc_GeVfm;
      s_momentum(i, 1) = py/hbarc_GeVfm;
      s_momentum(i, 2) = pz/hbarc_GeVfm;
      s_energy(i) = E/hbarc_GeVfm;
      s_mass(i) = mass;

      // discard charges
      s_baryon_charge(i) = 0.0;
      s_strangeness_charge(i) = 0.0;
      s_electric_charge(i) = 0.0;

      s_time(i) = t;

      // PID not meaningful for M-lines (scattering centers); keep a sentinel
      s_PID(i) = 0;
    }

    std::cout << "Loaded " << n_sources << " JEWEL M-line sources (sorted by time)" << std::endl;
    //print all sources
    for (std::size_t i = 0; i < n_sources; ++i)
    {
      std::cout << "Source particle " << i << ": "
                << "Position: (" << s_position(i,0) << ", " << s_position(i,1) << ", " << s_position(i,2) << "), "
                << "Momentum: (" << s_momentum(i,0) << ", " << s_momentum(i,1) << ", " << s_momentum(i,2) << "), "
                << "Energy: " << s_energy(i) << ", "
                << "Mass: " << s_mass(i) << ", "
                << "Baryon Charge: " << s_baryon_charge(i) << ", "
                << "Strangeness Charge: " << s_strangeness_charge(i) << ", "
                << "Electric Charge: " << s_electric_charge(i) << ", "
                << "Time: " << s_time(i)  << ", "
                << "PID: "  << s_PID(i)  << std::endl;
    }
  }


};
}// namespace ccake



#endif // SOURCE_H
