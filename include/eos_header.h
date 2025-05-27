#ifndef EOS_HEADER_H
#define EOS_HEADER_H

static constexpr bool use_rootfinder                      = true;
static constexpr bool use_static_C_library                = false;
static constexpr bool use_nonconformal_extension          = false;
static constexpr bool use_tanh_conformal                  = true;
static constexpr bool skip_failed_EoS                     = true;
static constexpr bool zero_unsolvable_charge_densities    = false;
static constexpr bool prohibit_unstable_cs2               = true;
static constexpr bool prohibit_acausal_cs2                = true;
static constexpr bool restrict_mu_T_ratios                = true;


static constexpr int VERBOSE = 3;
static constexpr double TINY = 1e-25;
static constexpr double TBQS_INFINITY = 1e10;  // indefinitely large limits for EoS

#endif
