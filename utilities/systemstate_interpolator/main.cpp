#include "system_state_int.h"
#include <iostream>
#include <Kokkos_Core.hpp>
#include <yaml-cpp/yaml.h>

int main(int argc, char** argv)
{
    Kokkos::initialize(argc, argv);
    int ret = 0;
    try {
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " <config.yaml>\n";
            ret = 1;
        } else {
            const YAML::Node cfg = YAML::LoadFile(argv[1]);
            const YAML::Node& ss = cfg["systemstate"];
            const YAML::Node& qc = cfg["query"];

            std::string output_dir = ss["output_dir"].as<std::string>();
            int    D            = ss["D"].as<int>();
            double t0           = ss["t0"].as<double>();
            double dt           = ss["dt"].as<double>();
            int    output_steps = ss["output_steps"].as<int>();
            double h_cell       = ss["h_cell"].as<double>(0.5);

            double t_query      = qc["t"].as<double>();
            double hT           = qc["hT"].as<double>();
            std::vector<double> xv = qc["x"].as<std::vector<double>>();
            double x[3] = {xv[0], xv[1], xv.size() > 2 ? xv[2] : 0.0};

            SystemState state(D, output_dir, t0, dt, output_steps);

            state.load_all(h_cell);
            std::cout << "Loaded " << state.evolution.size() << " snapshots\n";
            std::cout << "  time range: ["
                      << state.evolution.front().time << ", "
                      << state.evolution.back().time  << "]\n";

            auto val = state.interpolate_at(t_query, x, hT);

            std::cout << "Interpolated values at t=" << t_query
                      << " (" << x[0] << ", " << x[1] << ", " << x[2] << "):\n"
                      << "  sigma_star = " << val.sigma_star << "\n"
                      << "  e     = " << val.e     << "\n"
                      << "  s     = " << val.s     << "\n"
                      << "  rhoB  = " << val.rhoB  << "\n"
                      << "  rhoS  = " << val.rhoS  << "\n"
                      << "  rhoQ  = " << val.rhoQ  << "\n"
                      << "  T     = " << val.T     << "\n";
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        ret = 2;
    }
    Kokkos::finalize();
    return ret;
}
