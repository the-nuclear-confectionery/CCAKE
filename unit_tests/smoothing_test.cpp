#include <gtest/gtest.h>
#include <memory>
#include <Kokkos_Core.hpp>
//#include <gmock/gmock.h> //To mock settings pointers

#include "sph_workstation.h"
#include "system_state.h"
#include "settings.h"
#include "particle.h"
#define TOL 1.E-14

//class MockSettings : Settings{
//    MOCK_METHOD(void, check_consistency, (), (override));
//};

TEST(SmoothingTest, Smooth2D){

    int argc = 1;
    char* argv[] = {"smooth_test"};
    Kokkos::ScopeGuard guard(argc, argv);
    //Create mocked settings
    std::shared_ptr<Settings> settingsPtr = std::make_shared<Settings> ();
    settingsPtr->xmin = -4;
    settingsPtr->ymin = -6;
    settingsPtr->etamin = -2;
    settingsPtr->hT = 1;
    std::shared_ptr<ccake::SystemState<2>> sysPtr = std::make_shared<ccake::SystemState<2>> (settingsPtr);
    ccake::SPHWorkstation<2,ccake::EoM_default> ws(settingsPtr, sysPtr);

    //ccake::SystemState<2> sys(settingsPtr);
    sysPtr->hT = settingsPtr->hT;

    //Place particles in positions {0,0}, {.5,0}, {.8,.6} and {0,1.5}
    std::array<ccake::Particle<2>,4> particles_array_in;
    particles_array_in[0].r = {0,0};
    particles_array_in[1].r = {0,.5};
    particles_array_in[2].r = {.8,.6};
    particles_array_in[3].r = {1.5,0};
    for (auto p : particles_array_in) sysPtr->add_particle(p);
    sysPtr->n_particles = 4;
    for (int i=0; i<4; ++i){
        sysPtr->particles[i].ID = i;
        sysPtr->particles[i].sph_mass.s = 1;
        sysPtr->particles[i].sph_mass.rhoB = 2;
        sysPtr->particles[i].sph_mass.rhoQ = 3;
        sysPtr->particles[i].sph_mass.rhoS = 4;

        sysPtr->particles[i].extensive.s = 1;
        sysPtr->particles[i].extensive.rhoB = 2;
        sysPtr->particles[i].extensive.rhoQ = 3;
        sysPtr->particles[i].extensive.rhoS = 4;
    }
    
    sysPtr->allocate_cabana_particles();
    sysPtr->reset_neighbour_list();    
    ws.smooth_all_particle_fields(1.);
    sysPtr->copy_device_to_host();

    EXPECT_LT( fabs(sysPtr->particles[0].smoothed.s - 20./7./M_PI), TOL) << "Test failed! Reporting values used: " << endl
    << "Particle 1: -------" << endl
    << "r: " << sysPtr->particles[0].r << endl
    << "sph_mass: " << sysPtr->particles[0].sph_mass.s << endl
    << "extensive: " << sysPtr->particles[0].extensive.s << endl
    << "smoothed: " << sysPtr->particles[0].smoothed.s << endl
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

    EXPECT_LT( fabs(sysPtr->particles[0].smoothed.rhoB - 4*20./7./M_PI), TOL) << "Test failed! Reporting values used: " << endl
    << "Particle 1: -------" << endl
    << "r: " << sysPtr->particles[0].r << endl
    << "sph_mass: " << sysPtr->particles[0].sph_mass.rhoB << endl
    << "extensive: " << sysPtr->particles[0].extensive.rhoB << endl
    << "smoothed: " << sysPtr->particles[0].smoothed.rhoB << endl
    << "Particle 2: -------" << endl
    << "r: " << sysPtr->particles[1].r << endl
    << "sph_mass: " << sysPtr->particles[1].sph_mass.rhoB << endl
    << "extensive: " << sysPtr->particles[1].extensive.rhoB << endl
    << "smoothed: " << sysPtr->particles[1].smoothed.rhoB << endl
    << "Particle 3: -------" << endl
    << "r: " << sysPtr->particles[2].r << endl
    << "sph_mass: " << sysPtr->particles[2].sph_mass.rhoB << endl
    << "extensive: " << sysPtr->particles[2].extensive.rhoB << endl
    << "smoothed: " << sysPtr->particles[2].smoothed.rhoB << endl
    << "Particle 4: -------" << endl
    << "r: " << sysPtr->particles[3].r << endl
    << "sph_mass: " << sysPtr->particles[3].sph_mass.rhoB << endl
    << "extensive: " << sysPtr->particles[3].extensive.rhoB << endl
    << "smoothed: " << sysPtr->particles[3].smoothed.rhoB << endl;

    EXPECT_LT( fabs(sysPtr->particles[0].smoothed.rhoQ - 9*20./7./M_PI), TOL) << "Test failed! Reporting values used: " << endl
    << "Particle 1: -------" << endl
    << "r: " << sysPtr->particles[0].r << endl
    << "sph_mass: " << sysPtr->particles[0].sph_mass.rhoQ << endl
    << "extensive: " << sysPtr->particles[0].extensive.rhoQ << endl
    << "smoothed: " << sysPtr->particles[0].smoothed.rhoQ << endl
    << "Particle 2: -------" << endl
    << "r: " << sysPtr->particles[1].r << endl
    << "sph_mass: " << sysPtr->particles[1].sph_mass.rhoQ << endl
    << "extensive: " << sysPtr->particles[1].extensive.rhoQ << endl
    << "smoothed: " << sysPtr->particles[1].smoothed.rhoQ << endl
    << "Particle 3: -------" << endl
    << "r: " << sysPtr->particles[2].r << endl
    << "sph_mass: " << sysPtr->particles[2].sph_mass.rhoQ << endl
    << "extensive: " << sysPtr->particles[2].extensive.rhoQ << endl
    << "smoothed: " << sysPtr->particles[2].smoothed.rhoQ << endl
    << "Particle 4: -------" << endl
    << "r: " << sysPtr->particles[3].r << endl
    << "sph_mass: " << sysPtr->particles[3].sph_mass.rhoQ << endl
    << "extensive: " << sysPtr->particles[3].extensive.rhoQ << endl
    << "smoothed: " << sysPtr->particles[3].smoothed.rhoQ << endl;

    EXPECT_LT( fabs(sysPtr->particles[0].smoothed.rhoS - 16*20./7./M_PI), TOL) << "Test failed! Reporting values used: " << endl
    << "Particle 1: -------" << endl
    << "r: " << sysPtr->particles[0].r << endl
    << "sph_mass: " << sysPtr->particles[0].sph_mass.rhoS << endl
    << "extensive: " << sysPtr->particles[0].extensive.rhoS << endl
    << "smoothed: " << sysPtr->particles[0].smoothed.rhoS << endl
    << "Particle 2: -------" << endl
    << "r: " << sysPtr->particles[1].r << endl
    << "sph_mass: " << sysPtr->particles[1].sph_mass.rhoS << endl
    << "extensive: " << sysPtr->particles[1].extensive.rhoS << endl
    << "smoothed: " << sysPtr->particles[1].smoothed.rhoS << endl
    << "Particle 3: -------" << endl
    << "r: " << sysPtr->particles[2].r << endl
    << "sph_mass: " << sysPtr->particles[2].sph_mass.rhoS << endl
    << "extensive: " << sysPtr->particles[2].extensive.rhoS << endl
    << "smoothed: " << sysPtr->particles[2].smoothed.rhoS << endl
    << "Particle 4: -------" << endl
    << "r: " << sysPtr->particles[3].r << endl
    << "sph_mass: " << sysPtr->particles[3].sph_mass.rhoS << endl
    << "extensive: " << sysPtr->particles[3].extensive.rhoS << endl
    << "smoothed: " << sysPtr->particles[3].smoothed.rhoS << endl;

    EXPECT_LT( fabs(sysPtr->particles[0].hydro.sigma - 20./7./M_PI), TOL) << "Test failed! Reporting values used: " << endl
    << "Particle 1: -------" << endl
    << "r: " << sysPtr->particles[0].r << endl
    << "sph_mass: " << sysPtr->particles[0].sph_mass.s << endl
    << "extensive: " << sysPtr->particles[0].extensive.s << endl
    << "smoothed: " << sysPtr->particles[0].smoothed.s << endl
    << "sigma: " << sysPtr->particles[0].hydro.sigma << endl
    << "Particle 2: -------" << endl
    << "r: " << sysPtr->particles[1].r << endl
    << "sph_mass: " << sysPtr->particles[1].sph_mass.s << endl
    << "extensive: " << sysPtr->particles[1].extensive.s << endl
    << "smoothed: " << sysPtr->particles[1].smoothed.s << endl
    << "sigma: " << sysPtr->particles[1].hydro.sigma << endl
    << "Particle 3: -------" << endl
    << "r: " << sysPtr->particles[2].r << endl
    << "sph_mass: " << sysPtr->particles[2].sph_mass.s << endl
    << "extensive: " << sysPtr->particles[2].extensive.s << endl
    << "smoothed: " << sysPtr->particles[2].smoothed.s << endl
    << "sigma: " << sysPtr->particles[2].hydro.sigma << endl
    << "Particle 4: -------" << endl
    << "r: " << sysPtr->particles[3].r << endl
    << "sph_mass: " << sysPtr->particles[3].sph_mass.s << endl
    << "extensive: " << sysPtr->particles[3].extensive.s << endl
    << "smoothed: " << sysPtr->particles[3].smoothed.s << endl
    << "sigma: " << sysPtr->particles[3].hydro.sigma << endl;

}
