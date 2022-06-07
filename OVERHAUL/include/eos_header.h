#ifndef EOS_HEADER_H
#define EOS_HEADER_H

static constexpr bool use_rootfinder                      = true;
static constexpr bool use_delaunay                        = !use_rootfinder;
static constexpr bool use_static_C_library                = true;
//static constexpr bool accept_nearest_neighbor             = false;
//static constexpr bool discard_unsolvable_charge_densities = false;
//static constexpr bool use_conformal_as_fallback           = true;
static constexpr bool use_nonconformal_extension          = false;
static constexpr bool use_tanh_conformal                  = true;
static constexpr bool skip_failed_EoS                     = true;
static constexpr bool zero_unsolvable_charge_densities    = false;

static constexpr int VERBOSE = 3;
static constexpr double TINY = 1e-25;
static constexpr double TBQS_INFINITY = 1e10;  // indefinitely large limits for EoS

#endif