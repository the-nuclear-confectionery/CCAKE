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
    if (settingsPtr->source_type != "analytic"){
      //read source file
      std::cout << "Reading source file" << std::endl;
      read_source_file();
    }
    else{
      //initialize source array
      n_sources = 1.0;
      std::cout << "Initializing analytic source term" << std::endl;
      sources = source_array("sources", n_sources);
    }

  }


  void read_source_file()
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

      // Transform if hyperbolic and z² > t²
      if (settingsPtr->coordinate_system == "hyperbolic")
      {
        double tau = std::sqrt(t*t - z*z);
        double eta_s = 0.5 * std::log((t + z) / (t - z));

        double mT = std::sqrt(px*px + py*py + mass*mass);
        double y_particle = 0.5 * std::log((E + pz) / (E - pz));

        double p_tau = mT * std::cosh(y_particle - eta_s);
        double p_eta = mT * std::sinh(y_particle - eta_s);
        std::cout << "Transforming to hyperbolic coordinates" << std::endl;
        // Overwrite momentum & coordinates
        pz = p_eta;
        E  = p_tau;
        t  = tau;
        z  = eta_s;
      }
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


  void add_source()
  {
    CREATE_VIEW(device_, systemPtr->cabana_particles);
    auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

    double hT = settingsPtr->hT;
    double t = systemPtr->t;
    double dt = settingsPtr->dt;
    double t_prev = t - dt;
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
        double normalization = 1.; // IC normalization factor
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
      }
    }

    Kokkos::fence();
    std::cout << "Added source terms for " << n_sources_added << " source particles out of " << n_sources << std::endl;
  }
};
}// namespace ccake



#endif // SOURCE_H