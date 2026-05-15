#ifndef EOS_HEADER_H
#define EOS_HEADER_H

static constexpr bool use_rootfinder                      = true;
static constexpr bool use_static_C_library                = false;
static constexpr bool use_nonconformal_extension          = false;
static constexpr bool use_tanh_conformal                  = true;
static constexpr bool skip_failed_EoS                     = true;
static constexpr bool zero_unsolvable_charge_densities    = false;
// prohibit_unstable_cs2 and prohibit_acausal_cs2 are now runtime settings
// in Settings::prohibit_unstable_cs2 / prohibit_acausal_cs2 (see settings.h)
static constexpr bool restrict_mu_T_ratios                = false;


static constexpr int VERBOSE = 0;
static constexpr double TINY = 1e-25;
static constexpr double TBQS_INFINITY = 1e10;  // indefinitely large limits for EoS

#endif