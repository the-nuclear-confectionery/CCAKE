#ifndef EOS_HEADER_H
#define EOS_HEADER_H

static constexpr bool use_rootfinder                      = true;
static constexpr bool use_delaunay                        = !use_rootfinder;
static constexpr bool use_static_C_library                = true;
//static constexpr bool accept_nearest_neighbor             = false;
//static constexpr bool discard_unsolvable_charge_densities = false;
//static constexpr bool use_conformal_as_fallback           = true;
static constexpr bool use_nonconformal_extension          = true;

static constexpr int VERBOSE = 0;
static constexpr double TINY = 1e-25;
static constexpr double TBQS_INFINITY = 1e10;  // indefinitely large limits for EoS

#endif