//============================================================================//
// Header file for all mathematical and physical constants
//============================================================================//
#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace constants
{
  constexpr double pi          = 3.1415926535897932384626433;
  constexpr double twopi       = 2.0*pi;

  //constexpr double hbarc_GeVfm = 0.1973269804593;
  //constexpr double hbarc_MeVfm = 197.3269804593;
  constexpr double hbarc_GeVfm = 0.1973;
  constexpr double hbarc_MeVfm = 197.3;
  constexpr double hbarc       = hbarc_GeVfm;    // default to GeV units
  constexpr double hbarc2      = hbarc*hbarc;
  constexpr double hbarc3      = hbarc*hbarc*hbarc;
  constexpr double hbarcm1     = 1.0/hbarc;
  constexpr double hbarcm2     = 1.0/hbarc2;
  constexpr double hbarcm3     = 1.0/hbarc3;
	
  constexpr double MeVtoGeV    = 0.001;
  constexpr double GeVtoMeV    = 1000.0;
  constexpr double MeVtofm     = hbarc_MeVfm;
  constexpr double fmtoMeV     = 1.0/hbarc_MeVfm;
  constexpr double GeVtofm     = hbarc_GeVfm;
  constexpr double fmtoGeV     = 1.0/hbarc_GeVfm;
}

#endif