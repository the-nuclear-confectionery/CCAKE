#include <gtest/gtest.h>
#include <memory>
#include <Kokkos_Core.hpp>
//#include <gmock/gmock.h> //To mock settings pointers

#include "sph_workstation.h"
#include "system_state.h"
#include "settings.h"
#include "particle.h"
#include "eom_default.h"
#define TOL 5.E-13


std::shared_ptr<ccake::SystemState<2>> prepare_system(){
    std::shared_ptr<Settings> settingsPtr = std::make_shared<Settings> ();
    settingsPtr->xmin = -4;
    settingsPtr->ymin = -6;
    settingsPtr->etamin = -2;
    settingsPtr->hT = 1;
    settingsPtr->using_shear = true;

    std::shared_ptr<ccake::SystemState<2>> sysPtr = std::make_shared<ccake::SystemState<2>> (settingsPtr);
    sysPtr->hT = settingsPtr->hT;
    const int n_sph_particles = 1;
    std::array<ccake::Particle<2>,n_sph_particles> particles_array_in;
    particles_array_in[0].r = {-1.5,2.5};
    for (auto p : particles_array_in) sysPtr->add_particle(p);
    sysPtr->n_particles = n_sph_particles;
    sysPtr->t = 2.;

    return sysPtr;
}

TEST(eom_default,get_LRF){
    
    ASSERT_DOUBLE_EQ(ccake::EoM_default<2>::get_LRF(4.,2.,2.),1.);
}

TEST(eom_default, dsigma_dt){
    

    std::shared_ptr<ccake::SystemState<2>> systemPtr = prepare_system();
    //Setup test
    systemPtr->particles[0].ID = 0;
    systemPtr->particles[0].hydro.sigma = 2.;
    systemPtr->particles[0].hydro.gradV(0,0) = 1.;
    systemPtr->particles[0].hydro.gradV(1,1) = 2.;
    systemPtr->particles[0].hydro.gradV(0,1) = 4.5; //Should not matter
    systemPtr->particles[0].hydro.gradV(1,0) = 4.5; //Should not matter
    double expected_dsigma_dt = -6;
    systemPtr->allocate_cabana_particles();
    ccake::EoM_default<2>::evaluate_time_derivatives(systemPtr);
    systemPtr->copy_device_to_host();
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.dsigma_dt, expected_dsigma_dt);
}

TEST(eom_default, gamma_and_v){
    std::shared_ptr<ccake::SystemState<2>> systemPtr = prepare_system();
    //Setup test
    systemPtr->particles[0].ID = 0;
    systemPtr->particles[0].hydro.sigma = 2.;
    systemPtr->particles[0].hydro.u(0) = 3.;
    systemPtr->particles[0].hydro.u(1) = 4.;
    double expected_gamma = sqrt(26.0);
    double expected_vx = 3./sqrt(26.0);
    double expected_vy = 4./sqrt(26.0);

    systemPtr->particles[0].hydro.sigma = 2.;
    systemPtr->particles[0].hydro.gradV(0,0) = 1.;
    systemPtr->particles[0].hydro.gradV(1,1) = 2.;
    systemPtr->particles[0].hydro.gradV(0,1) = 4.5; //Should not matter
    systemPtr->particles[0].hydro.gradV(1,0) = 4.5; //Should not matter
    double expected_dsigma_dt = -2.0 * (1.0 + 2.0);

    systemPtr->allocate_cabana_particles();
    ccake::EoM_default<2>::reset_pi_tensor(systemPtr);
    systemPtr->copy_device_to_host();
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.gamma, expected_gamma);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.v(0), expected_vx);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.v(1), expected_vy);
}

TEST(eom_default, test_time_derivative){
    std::shared_ptr<ccake::SystemState<2>> systemPtr = prepare_system();
    //Setup test
    
    systemPtr->particles[0].ID = 0.;
    systemPtr->particles[0].hydro.sigma = 2.;
    systemPtr->particles[0].hydro.u(0) = 3.;
    systemPtr->particles[0].hydro.u(1) = 4.;
    systemPtr->particles[0].hydro.Bulk = 6.;
    systemPtr->particles[0].hydro.setas = 2.*M_PI;
    systemPtr->particles[0].hydro.stauRelax = 3.*M_PI;
    systemPtr->particles[0].hydro.zeta = 4.*M_PI;
    systemPtr->particles[0].hydro.tauRelax = 2.*M_PI;
    systemPtr->particles[0].thermo.s = 7.;
    systemPtr->particles[0].thermo.T = M_PI/2.;
    systemPtr->particles[0].thermo.w = M_PI/2.;
    systemPtr->particles[0].thermo.dwds = M_PI;
    systemPtr->particles[0].thermo.dwdB = 8.;
    systemPtr->particles[0].thermo.dwdQ = 9.;
    systemPtr->particles[0].thermo.dwdS = 10.;
    systemPtr->particles[0].thermo.rhoB = 11.;
    systemPtr->particles[0].thermo.rhoQ = 12.;
    systemPtr->particles[0].thermo.rhoS = 13.;
    systemPtr->particles[0].hydro.gradshear = {7./2.,9./2.};
    systemPtr->particles[0].hydro.gradBulk = {7./5.,9./5.};
    systemPtr->particles[0].hydro.gradP = {7./4.,9./4.};
    systemPtr->particles[0].hydro.divshear = {7./8.,9./8.};

    systemPtr->particles[0].hydro.sigma = 2.;
    systemPtr->particles[0].hydro.gradV(0,0) = 1.;
    systemPtr->particles[0].hydro.gradV(1,1) = 2.;
    systemPtr->particles[0].hydro.gradV(0,1) = 4.5; //Should not matter
    systemPtr->particles[0].hydro.gradV(1,0) = 5.5; //Should not matter

    systemPtr->particles[0].hydro.shv(1,1) = 3.;
    systemPtr->particles[0].hydro.shv(1,2) = 4.;
    systemPtr->particles[0].hydro.shv(2,2) = 5.;
    //systemPtr->particles[0].hydro.shv(2,1) = 4.;
    //systemPtr->particles[0].hydro.shv(0,1) = 25./sqrt(26.);
    //systemPtr->particles[0].hydro.shv(1,0) = 25./sqrt(26.);
    //systemPtr->particles[0].hydro.shv(0,2) = 32./sqrt(26.);
    //systemPtr->particles[0].hydro.shv(2,0) = 32./sqrt(26.);
    //systemPtr->particles[0].hydro.shv(0,0) = 203./26;
    //systemPtr->particles[0].hydro.shv33 = -5./104.;
    //systemPtr->particles[0].hydro.gamma = sqrt(26.);
    //systemPtr->particles[0].hydro.gamma_squared = 26.;
    //systemPtr->particles[0].hydro.gamma_cube = 26.*sqrt(26.);
    //systemPtr->particles[0].hydro.t = 2;
    //systemPtr->particles[0].hydro.gamma_tau = sqrt(26.)*2;
    //systemPtr->particles[0].hydro.v(0) = systemPtr->particles[0].hydro.u(0)/sqrt(26.);
    //systemPtr->particles[0].hydro.v(1) = systemPtr->particles[0].hydro.u(1)/sqrt(26.);

    double expected_gamma = sqrt(26.0);
    double expected_gamma_squared = 26.0;
    double expected_gamma_cube = 26.0*sqrt(26.0);
    double expected_gamma_tau = 2.*sqrt(26.0);
    double expected_vx = 3./sqrt(26.0);
    double expected_vy = 4./sqrt(26.0);
    double expected_dwdsT1 = -1;
    double expected_dsigma_dt = -6;
    double expected_sigl = -3.5;
    double expected_bulk = 3.*sqrt(2/13.0);
    double expected_C = M_PI/2.+3.*sqrt(2/13.0);
    double expected_eta_o_tau = 2./3.;
    double expected_Agam = -328.-6.*sqrt(2./13.)-13.*M_PI/2.;
    double expected_C_tot = M_PI/2.+3.*sqrt(2/13.0)-25./39;
    double expected_Agam2 = -6./13. -74977./234./sqrt(26.) - sqrt(13./2.)*M_PI/2.;
    double expected_bsub = 69./sqrt(26.);
    double expected_Btot = 42. + 269465./9./sqrt(26.)+3./sqrt(26.)/M_PI+91.*M_PI*sqrt(13./2.)/2;

    double expected_du_dt0=-9.1234319302411276;
    double expected_du_dt1=-16.259788577386033;
    double expected_inside = -15.534563014702962;
    double expected_d_dt_extensive_s = -4.7377394923567692;
    double expected_dBulk_dt = 0.36555469870526137;
    double expected_dshv_dt00 = -86.251382196576950;
    double expected_dshv_dt01 = -85.958565490964447;
    double expected_dshv_dt10 = -85.958565490964447;
    double expected_dshv_dt11 = -63.230392427359653;


    double expected_uu00 = 9.;
    double expected_uu01 = 12.;
    double expected_uu10 = 12.;
    double expected_uu11 = 16.;
    double expected_pimin00 = 3.;
    double expected_pimin01 = 4.;
    double expected_pimin10 = 4.;
    double expected_pimin11 = 5.;
    double expected_piu00 = 75./sqrt(26.);
    double expected_piu01 = 48.*sqrt(2./13.);
    double expected_piu10 = 50.*sqrt(2./13.);
    double expected_piu11 = 64.*sqrt(2./13.);
    double expected_piutot00 = 75.*sqrt(2./13.);
    double expected_piutot01 = 98.*sqrt(2./13.);
    double expected_piutot10 = 98.*sqrt(2./13.);
    double expected_piutot11 = 128.*sqrt(2./13.);
    double expected_div_u = -2.825925663130984;
    double expected_extensivetheta = -0.55283181266918218;
    
    
    
    double expected_shv01 = 25/expected_gamma;
    double expected_shv02 = 32/expected_gamma;
    double expected_shv00 = 203./26.;
    double expected_shv33 =-5./104.;

    double expected_gradU00 = 76.*sqrt(26.);
    double expected_gradU01 = 69.*sqrt(26.);
    double expected_gradU10 = 211.*sqrt(13./2.);
    double expected_gradU11 = 88.*sqrt(26.);
    

    //Execute the code
    systemPtr->allocate_cabana_particles();
    ccake::EoM_default<2>::reset_pi_tensor(systemPtr);
    ccake::EoM_default<2>::evaluate_time_derivatives(systemPtr);
    systemPtr->copy_device_to_host();


    //Asserts that we have the correct shear tensor
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.shv(2,1), systemPtr->particles[0].hydro.shv(1,2));
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.shv(0,1), expected_shv01);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.shv(1,0), expected_shv01);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.shv(0,2), expected_shv02);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.shv(2,0), expected_shv02);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.shv(0,0), expected_shv00);
    ASSERT_LE(std::fabs(systemPtr->particles[0].hydro.shv33-expected_shv33), TOL);

    //Assert that we correctly computed the evolution variables
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.dsigma_dt, expected_dsigma_dt);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.gamma, expected_gamma);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.v(0), expected_vx);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.v(1), expected_vy);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.gamma_squared, expected_gamma_squared);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.gamma_cube, expected_gamma_cube);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.gamma_tau, expected_gamma_tau);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.dwdsT1, expected_dwdsT1);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.dsigma_dt, expected_dsigma_dt);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.sigl, expected_sigl);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.extensivePI, expected_bulk);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.C, expected_C);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.eta_o_tau, expected_eta_o_tau);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.Agam, expected_Agam);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.Agam2, expected_Agam2);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.Btot, expected_Btot);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.Btot, expected_Btot);


    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.uu(0,0), expected_uu00);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.uu(0,1), expected_uu10);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.uu(1,0), expected_uu01);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.uu(1,1), expected_uu11);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.piu(0,0), expected_piu00);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.piu(0,1), expected_piu10);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.piu(1,0), expected_piu01);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.piu(1,1), expected_piu11);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.piutot(0,0), expected_piutot00);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.piutot(0,1), expected_piutot10);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.piutot(1,0), expected_piutot01);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.piutot(1,1), expected_piutot11);
    
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.gradU(0,0), expected_gradU00);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.gradU(1,1), expected_gradU11);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.gradU(0,1), expected_gradU01);
    ASSERT_DOUBLE_EQ(systemPtr->particles[0].hydro.gradU(1,0), expected_gradU10);
    
    ASSERT_LE( fabs(systemPtr->particles[0].hydro.du_dt(0)/expected_du_dt0 - 1), TOL) ;
    ASSERT_LE( fabs(systemPtr->particles[0].hydro.du_dt(1)/expected_du_dt1 - 1), TOL) ;
    
    ASSERT_LE( fabs(systemPtr->particles[0].hydro.div_u/expected_div_u - 1), TOL) ;
    ASSERT_LE( fabs(systemPtr->particles[0].hydro.extensivetheta/expected_extensivetheta - 1), TOL) ;
    ASSERT_LE( fabs(systemPtr->particles[0].hydro.inside/expected_inside - 1), TOL) ;

    ASSERT_LE( fabs(systemPtr->particles[0].d_dt_extensive.s/expected_d_dt_extensive_s - 1), TOL) ;
    ASSERT_LE( fabs(systemPtr->particles[0].hydro.dBulk_dt/expected_dBulk_dt - 1), TOL) ;

    ASSERT_LE( fabs(systemPtr->particles[0].hydro.dshv_dt(0,0)/expected_dshv_dt00 - 1), TOL) ;
    ASSERT_LE( fabs(systemPtr->particles[0].hydro.dshv_dt(1,0)/expected_dshv_dt01 - 1), TOL) ;
    ASSERT_LE( fabs(systemPtr->particles[0].hydro.dshv_dt(0,1)/expected_dshv_dt10 - 1), TOL) ;
    ASSERT_LE( fabs(systemPtr->particles[0].hydro.dshv_dt(1,1)/expected_dshv_dt11 - 1), TOL) ;



}

// Define the main function
int main(int argc, char* argv[]) {
 
    // Initialize Kokkos
    Kokkos::ScopeGuard guard(argc, argv);

    // Run the tests
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}