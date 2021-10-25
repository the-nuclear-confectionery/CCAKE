#include <cmath>
#include "vector.h"
#include "matrix.h"
#include "LinkList.h"
#include "particle.h"
#include "eos.h"
#include "runge_kutta.h"

namespace RK
{
  void bsq_second_order( double dx, void (*derivatives)( LinkList &linklist ),
                         LinkList &linklist )
  {
    int N = linklist.n();

    // creating vectors of vectors of the derivatives at each step
    vector<double> etasigma0(N);
    vector<double> Bulk0(N);

    vector< Vector<double,2> > u0(N);
    vector< Vector<double,2> > r0(N);

    vector< Matrix <double,2,2> > shv0(N);

    double E0, t0;

    linklist.rk2 = 1;
    t0           = linklist.t;

    // initialize quantities at current time step
    for (int i=0; i<N; ++i)
    {
      u0[i]        = linklist._p[i].u;
      r0[i]        = linklist._p[i].r;
      etasigma0[i] = linklist._p[i].eta_sigma;
      Bulk0[i]     = linklist._p[i].Bulk;
      mini( shv0[i], linklist._p[i].shv );
    }

    E0 = linklist.Ez;

    ////////////////////////////////////////////
    //    first step
    ////////////////////////////////////////////

    // compute derivatives
    (*derivatives)(linklist);

    // update quantities
    for (int i=0; i<N; ++i)
    {
      const auto & this_particle = linklist._p[i];
      this_particle.u            = u0[i]        + 0.5*dx*this_particle.du_dt;
      this_particle.r            = r0[i]        + 0.5*dx*this_particle.v;
      this_particle.eta_sigma    = etasigma0[i] + 0.5*dx*this_particle.detasigma_dt;
      this_particle.Bulk         = Bulk0[i]     + 0.5*dx*this_particle.dBulk_dt;
      tmini( this_particle.shv,    shv0[i]      + 0.5*dx*this_particle.dshv_dt );
    }

    linklist.Ez = E0 + 0.5*dx*linklist.dEz;
    linklist.t  = t0 + 0.5*dx;

    ////////////////////////////////////////////
    //    second step
    ////////////////////////////////////////////

    // compute derivatives
    (*derivatives)(linklist);

    // update quantities
    for (int i=0; i<N; ++i)
    {
      const auto & this_particle = linklist._p[i];
      this_particle.u            = u0[i]        + dx*this_particle.du_dt;
      this_particle.r            = r0[i]        + dx*this_particle.v;
      this_particle.eta_sigma    = etasigma0[i] + dx*this_particle.detasigma_dt;
      tmini( this_particle.shv,    shv0[i]      + dx*this_particle.dshv_dt );
    }

    linklist.Ez = E0 + dx*linklist.dEz;
    linklist.t  = t0 + dx;

	}

}