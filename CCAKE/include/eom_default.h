#ifndef EOM_DEFAULT_H
#define EOM_DEFAULT_H

#include "particle.h"
#include "utilities.h"
#include "densities.h"
#include "eom.h"
#include "hydrodynamic_info.h"
#include "settings.h"
#include "thermodynamic_info.h"
#include "vector.h"
#include "milne.hpp"
#include "system_state.h"

#include <Cabana_Core.hpp>

#include <memory>
namespace ccake{
template<unsigned int D>
class EoM_default
{
  private:
    static constexpr int VERBOSE = 3;
    std::shared_ptr<Settings> settingsPtr;
    Matrix<double,D,D> Imat;

  public:
    EoM_default() = delete;
    EoM_default( std::shared_ptr<Settings> settingsPtr_in ): settingsPtr(settingsPtr_in) { };
    virtual ~EoM_default(){};

    KOKKOS_FUNCTION
    static double gamma_calc(double u[D], const double &time_squared);
    KOKKOS_FUNCTION
    static double dot(double v[D], double u[D], const double &time_squared);
    KOKKOS_FUNCTION
    static void dot(double v[D], double T[D][D], const double &time_squared,
              double *result);
    KOKKOS_FUNCTION
    static double get_shvDD(double* pi_diag, const double &time_squared);
    KOKKOS_FUNCTION
    static double get_LRF(const double &lab, const double &gamma, const double &time_squared);
    KOKKOS_FUNCTION
    static double full_contraction(double M[D][D], double T[D][D], const double &time_squared);

    static void reset_pi_tensor(std::shared_ptr<SystemState<D>> sysPtr);
    static void evaluate_time_derivatives( std::shared_ptr<SystemState<D>> sysPtr);

    std::string name = "Israel-Stewart";                    // name associated to EoM
};
}
#endif