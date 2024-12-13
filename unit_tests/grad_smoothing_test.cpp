#include <gtest/gtest.h>
#include <memory>
#include <Kokkos_Core.hpp>
//#include <gmock/gmock.h> //To mock settings pointers

#include "sph_workstation.h"
#include "system_state.h"
#include "settings.h"
#include "particle.h"
#define TOL 1.E-14


TEST(GradSmoothingTest, Smooth2D){

    int argc = 1;
    char* argv[] = {"grad_smooth_test"};
    Kokkos::ScopeGuard guard(argc, argv);
    //Create mocked settings
    std::shared_ptr<Settings> settingsPtr = std::make_shared<Settings> ();
    settingsPtr->xmin = -4;
    settingsPtr->ymin = -6;
    settingsPtr->etamin = -2;
    settingsPtr->hT = 1;
    std::shared_ptr<ccake::SystemState<2>> sysPtr = std::make_shared<ccake::SystemState<2>> (settingsPtr);
    ccake::SPHWorkstation<2,ccake::EoM_default> ws(settingsPtr, sysPtr);
    const double target_val_x = 45./8./M_PI;
    const double target_val_y = 0;

    //ccake::SystemState<2> sys(settingsPtr);
    sysPtr->hT = settingsPtr->hT;

    //Place particles in positions {0,0}, {.5,0}, {.8,.6} and {0,1.5}
    const int n_sph_particles = 5;
    std::array<ccake::Particle<2>,n_sph_particles> particles_array_in;
    particles_array_in[0].r = {0,0}; //Should not contribute
    particles_array_in[1].r = {.5,0};
    particles_array_in[2].r = {1.0,0};
    particles_array_in[3].r = {1.5,0};
    particles_array_in[4].r = {-1.5,1.5}; //Should not contrinute
    //particles_array_in[4].r = {0,0};
    //particles_array_in[5].r = {-.5,0};
    //particles_array_in[6].r = {-1.0,0};
    //particles_array_in[7].r = {-1.5,0};
    for (auto p : particles_array_in) sysPtr->add_particle(p);
    sysPtr->n_particles = n_sph_particles;
    for (int i=0; i<n_sph_particles; ++i){
        sysPtr->particles[i].ID = i;
        sysPtr->particles[i].sph_mass.s = 1;
        sysPtr->particles[i].sph_mass.rhoB = 2;
        sysPtr->particles[i].sph_mass.rhoQ = 3;
        sysPtr->particles[i].sph_mass.rhoS = 4;
        sysPtr->particles[i].hydro.sigma = 1;
        sysPtr->particles[i].thermo.p = 1;


        sysPtr->particles[i].extensive.s = 1;
        sysPtr->particles[i].extensive.rhoB = 2;
        sysPtr->particles[i].extensive.rhoQ = 3;
        sysPtr->particles[i].extensive.rhoS = 4;
    }
    sysPtr->particles[1].thermo.p = .5;
    sysPtr->particles[2].thermo.p = 1.5;
    sysPtr->particles[3].thermo.p = 2.5;

    sysPtr->allocate_cabana_particles();
    sysPtr->reset_neighbour_list();    
    ws.smooth_all_particle_gradients(1.);
    sysPtr->copy_device_to_host();
    cout << "gradP is: " << sysPtr->particles[0].hydro.gradP(0) << endl;
    EXPECT_LT( fabs(sysPtr->particles[0].hydro.gradP(0) - target_val_x), TOL) << "Test failed! Reporting values used: " << endl
    << "Particle 1: -------" << endl
    << "r: " << sysPtr->particles[0].r << endl
    << "sph_mass: " << sysPtr->particles[0].sph_mass.s << endl
    << "extensive: " << sysPtr->particles[0].extensive.s << endl
    << "smoothed: " << sysPtr->particles[0].smoothed.s << endl
    << "sigma^*: " << sysPtr->particles[0].hydro.sigma << endl
    << "p: " << sysPtr->particles[0].thermo.p << endl
    << "grad_p: " << sysPtr->particles[0].hydro.gradP << endl
    << "Particle 2: -------" << endl
    << "r: " << sysPtr->particles[1].r << endl
    << "sph_mass: " << sysPtr->particles[1].sph_mass.s << endl
    << "extensive: " << sysPtr->particles[1].extensive.s << endl
    << "smoothed: " << sysPtr->particles[1].smoothed.s << endl
    << "Particle 3: -------" << endl
    << "r: " << sysPtr->particles[2].r << endl
    << "sph_mass: " << sysPtr->particles[2].sph_mass.s << endl
    << "extensive: " << sysPtr->particles[2].extensive.s << endl
    << "smoothed: " << sysPtr->particles[2].smoothed.s << endl
    << "Particle 4: -------" << endl
    << "r: " << sysPtr->particles[3].r << endl
    << "sph_mass: " << sysPtr->particles[3].sph_mass.s << endl
    << "extensive: " << sysPtr->particles[3].extensive.s << endl
    << "smoothed: " << sysPtr->particles[3].smoothed.s << endl;

 

}
