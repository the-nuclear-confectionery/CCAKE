#ifndef EoM_cartesian_H
#define EoM_cartesian_H

#include "particle.h"
#include "utilities.h"
#include "densities.h"
#include "hydrodynamic_info.h"
#include "settings.h"
#include "thermodynamic_info.h"
#include "vector.h"
#include "cartesian.hpp"
#include "system_state.h"

#include <Cabana_Core.hpp>

#include <memory>
namespace ccake{
/// @class EoM_cartesian
/// @brief  Equations of motion for the hydrodynamic evolution in Cartesian coordinates.
/// @details This class contains the default equations of motion 
/// (Israel-Stewart) for the hydrodynamic evolution.
/// @tparam D The number of spatial dimensions.
template<unsigned int D>
class EoM_cartesian
{
  private:
    std::shared_ptr<Settings> settingsPtr;

  public:
    EoM_cartesian() = delete;
    EoM_cartesian( std::shared_ptr<Settings> settingsPtr_in ): settingsPtr(settingsPtr_in) { };
    ~EoM_cartesian(){};

    KOKKOS_FUNCTION
    static double gamma_calc(double u[D], const double &time_squared);
    KOKKOS_FUNCTION
    static double get_LRF(const double &lab, const double &gamma, const double &time_squared);
    static void update_velocity(std::shared_ptr<SystemState<D>> sysPtr);
    static void reset_pi_tensor(std::shared_ptr<SystemState<D>> sysPtr);
    static void evaluate_time_derivatives( std::shared_ptr<SystemState<D>> sysPtr, std::shared_ptr<Settings> settingsPtr);
    static void calculate_MRF_shear(std::shared_ptr<SystemState<D>> sysPtr);
    static void compute_hydro_numbers(std::shared_ptr<SystemState<D>> sysPtr);
    std::string name = "Cartesian"; ///< name associated to EoM
};
}
#endif
