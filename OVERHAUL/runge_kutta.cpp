#include <cmath>
#include "vector.h"
#include "matrix.h"
#include "linklist.h"
#include "particle.h"
#include "eos.h"
#include "runge_kutta.h"

namespace RK
{
  void bsq_second_order( double dx, EquationsOfMotion & eom, SystemState & system )
  {
    int N = system.n();

    // creating vectors of vectors of the derivatives at each step
    vector<double> etasigma0(N);
    vector<double> Bulk0(N);

    vector< Vector<double,2> > u0(N);
    vector< Vector<double,2> > r0(N);

    vector< Matrix <double,2,2> > shv0(N);

    double E0, t0;

    system.rk2 = 1;
    t0           = system.t;

    // initialize quantities at current time step
    for (int i=0; i<N; ++i)
    {
      u0[i]        = system._p[i].u;
      r0[i]        = system._p[i].r;
      etasigma0[i] = system._p[i].eta_sigma;
      Bulk0[i]     = system._p[i].Bulk;
      mini( shv0[i], system._p[i].shv );
    }

    E0 = system.Ez;

    ////////////////////////////////////////////
    //    first step
    ////////////////////////////////////////////

    // compute derivatives
    //(*derivatives)(linklist);
    eom.BSQshear(system);

    // update quantities
    for (int i=0; i<N; ++i)
    {
      const auto & this_particle = system._p[i];
      this_particle.u            = u0[i]        + 0.5*dx*this_particle.du_dt;
      this_particle.r            = r0[i]        + 0.5*dx*this_particle.v;
      this_particle.eta_sigma    = etasigma0[i] + 0.5*dx*this_particle.detasigma_dt;
      this_particle.Bulk         = Bulk0[i]     + 0.5*dx*this_particle.dBulk_dt;
      tmini( this_particle.shv,    shv0[i]      + 0.5*dx*this_particle.dshv_dt );
    }

    system.Ez = E0 + 0.5*dx*system.dEz;
    system.t  = t0 + 0.5*dx;

    ////////////////////////////////////////////
    //    second step
    ////////////////////////////////////////////

    // compute derivatives
    //(*derivatives)(linklist);
    eom.BSQshear(system);

    // update quantities
    for (int i=0; i<N; ++i)
    {
      const auto & this_particle = system._p[i];
      this_particle.u            = u0[i]        + dx*this_particle.du_dt;
      this_particle.r            = r0[i]        + dx*this_particle.v;
      this_particle.eta_sigma    = etasigma0[i] + dx*this_particle.detasigma_dt;
      /*!!!!!!!!!!*/this_particle.Bulk         = Bulk0[i]     + 0.5*dx*this_particle.dBulk_dt;
      tmini( this_particle.shv,    shv0[i]      + dx*this_particle.dshv_dt );
    }

    system.Ez = E0 + dx*system.dEz;
    system.t  = t0 + dx;

	}

}