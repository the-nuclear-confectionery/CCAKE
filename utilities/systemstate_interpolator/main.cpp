#include "system_state_int.h"
#include <iostream>

int main(int argc, char** argv)
{
    try {
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " <output_directory>\n";
            return 1;
        }

        std::string output_dir = argv[1];

        // ---- hardcoded example values (or read from config) ----
        int    D            = 2;     // or 3
        double t0           = 0.6;
        double dt           = 0.05;
        int    output_steps = 1;

        // ---- construct and load system state ----
        SystemState state(
            D,
            output_dir,
            t0,
            dt,
            output_steps
        );

        double t = 1.;
        state.load(t);

        std::cout << "Loaded "
                  << state.particles.size()
                  << " particles at t = "
                  << state.file_time
                  << "\n";

        // ---- SPH interpolation at a point ----
        double x[3] = {-6.3, -1.48, 0.0};
        double hT   = 0.1;

        auto val = state.interpolate(x, hT);

        std::cout << "Interpolated values at (-6.3,-1.48,0):\n"
                  << "  sigma_star = " << val.sigma_star << "\n"
                  << "  e     = " << val.e << "\n"
                  << "  s     = " << val.s << "\n"
                  << "  rhoB  = " << val.rhoB << "\n"
                  << "  rhoS  = " << val.rhoS << "\n"
                  << "  rhoQ  = " << val.rhoQ << "\n"
                  << "  T     = " << val.T << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    }

    return 0;
}
