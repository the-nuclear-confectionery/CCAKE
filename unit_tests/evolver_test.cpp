#include <gtest/gtest.h>
#include <memory>
//#include <gmock/gmock.h> //To mock settings pointers

#include "sph_workstation.h"
#include "system_state.h"
#include "settings.h"
#include "particle.h"
#include "evolver.h"
#define TOL 5.E-13

namespace cc = ccake;

std::tuple<std::shared_ptr<cc::SystemState<2>>,std::shared_ptr<cc::Evolver<2>>> prepare_system(){
    std::shared_ptr<Settings> settingsPtr = std::make_shared<Settings> ();
    settingsPtr->xmin = -4;
    settingsPtr->ymin = -6;
    settingsPtr->etamin = -2;
    settingsPtr->hT = 1;
    settingsPtr->using_shear = true;
    settingsPtr->t0 = 2;

    std::shared_ptr<cc::SystemState<2>> sysPtr = std::make_shared<cc::SystemState<2>>(settingsPtr);
    sysPtr->n_particles=1;
    sysPtr->t=2.;
    cc::Particle<2> p;
    sysPtr->particles.push_back(p);
    sysPtr->initialize();

    std::shared_ptr<cc::Evolver<2>> evo = std::make_shared<cc::Evolver<2>>(settingsPtr,sysPtr); 

    return std::make_tuple(sysPtr,evo);
}

 
TEST(evolver, test_rk2){

    std::tuple<std::shared_ptr<cc::SystemState<2>>,std::shared_ptr<cc::Evolver<2>>> sys_evo_tuple;
    sys_evo_tuple = prepare_system();
    std::shared_ptr<cc::SystemState<2>> systemPtr = std::get<0>(sys_evo_tuple);
    std::shared_ptr<cc::Evolver<2>> evoPtr = std::get<1>(sys_evo_tuple);

    cc::Particle<2>* p = &systemPtr->particles[0];
    
    //Initial positions
    p->r(0) = -1.5;
    p->r(1) = -2.5;

    //Initial velocity
    p->hydro.u(0) = 3.;
    p->hydro.u(1) = 4.;
    p->hydro.gamma = sqrt(26.);
    p->hydro.v(0) = 3./sqrt(26.);
    p->hydro.v(1) = 4./sqrt(26.);

    //Initial entropy
    p->specific.s = 2.;
    
    //Initial Bulk
    p->hydro.Bulk = 6.;
    
    //Initial viscous shear tensor
    p->hydro.shv(1,1) = 3.;
    p->hydro.shv(1,2) = 4.;
    p->hydro.shv(2,2) = 5.;
    p->hydro.shv(2,1) = 4.;
    p->hydro.shv(0,1) = 25./sqrt(26.);
    p->hydro.shv(1,0) = 25./sqrt(26.);
    p->hydro.shv(0,2) = 32./sqrt(26.);
    p->hydro.shv(2,0) = 32./sqrt(26.);
    p->hydro.shv(0,0) = 203./26;
    p->hydro.shv33 = -5./104.;

    //Setup derivatives
    p->hydro.du_dt(0) = 9./2.;
    p->hydro.du_dt(1) = 7./2.;

    p->d_dt_spec.s = 5./2.;
    p->hydro.dBulk_dt = 3./2.;

    p->hydro.dshv_dt(0,0) = M_PI;
    p->hydro.dshv_dt(0,1) = 2*M_PI;
    p->hydro.dshv_dt(1,0) = 2*M_PI;
    p->hydro.dshv_dt(1,1) = 3*M_PI;

    p->ID = 0;
    
    double dt = .5;
    
    double const expected_r0 = -1.5 + dt*3./sqrt(26.);
    double const expected_r1 = -2.5 + dt*4./sqrt(26.);

    double const expected_u0 = 3 + dt*9./2.;
    double const expected_u1 = 4 + dt*7./2.;

    double const expected_specific_s = 2. + dt*5./2.;
    double const expected_Bulk = 6. + dt*3./2.;
    double const expected_shv11 = 3. + dt*M_PI;
    double const expected_shv12 = 4. + dt*2*M_PI;
    double const expected_shv21 = 4. + dt*2*M_PI;
    double const expected_shv22 = 5. + dt*3.*M_PI;

    auto evaluate_time_derivatives = [](){};

    systemPtr->allocate_cabana_particles();
    //systemPtr->copy_host_to_device();
    evoPtr->advance_timestep_rk2(dt, evaluate_time_derivatives);
    systemPtr->copy_device_to_host();
    
    //Asserts that we have the correct shear tensor
    ASSERT_DOUBLE_EQ(p->r(0), expected_r0);
    ASSERT_DOUBLE_EQ(p->r(1), expected_r1);
    ASSERT_DOUBLE_EQ(p->hydro.u(0), expected_u0);
    ASSERT_DOUBLE_EQ(p->hydro.u(1), expected_u1);
    ASSERT_DOUBLE_EQ(p->specific.s, expected_specific_s);
    ASSERT_DOUBLE_EQ(p->hydro.Bulk, expected_Bulk);
    
    ASSERT_DOUBLE_EQ(p->hydro.shv(1,1), expected_shv11);
    ASSERT_DOUBLE_EQ(p->hydro.shv(1,2), expected_shv12);
    ASSERT_DOUBLE_EQ(p->hydro.shv(2,1), expected_shv21);
    ASSERT_DOUBLE_EQ(p->hydro.shv(2,2), expected_shv22);
}

// Define the main function
int main(int argc, char* argv[]) {
 
    // Initialize Kokkos
    Kokkos::ScopeGuard guard(argc, argv);

    // Run the tests
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}