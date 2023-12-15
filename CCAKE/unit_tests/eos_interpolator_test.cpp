#include <gtest/gtest.h>
#include <memory>
#include <Cabana_Core.hpp>
//#include <gmock/gmock.h> //To mock settings pointers

#include "system_state.h"
#include "settings.h"
#include "particle.h"
#include "eos_interpolator.h"
#include "hydrodynamic_info.h"
#include "thermodynamic_info.h"
#include "densities.h"

#define TOL 5.E-13

namespace cc = ccake;

std::shared_ptr<cc::SystemState<2>> prepare_system(){
    int const num_particles = 10000;
    std::shared_ptr<Settings> settingsPtr = std::make_shared<Settings> ();
    settingsPtr->xmin = -4;
    settingsPtr->ymin = -6;
    settingsPtr->etamin = -2;
    settingsPtr->hT = 1;
    settingsPtr->using_shear = true;
    settingsPtr->t0 = 2;

    std::shared_ptr<cc::SystemState<2>> sysPtr = std::make_shared<cc::SystemState<2>>(settingsPtr);
    sysPtr->n_particles=num_particles;
    sysPtr->t=2.;
    for (int i=0; i<num_particles; i++){
        cc::Particle<2> p;
        sysPtr->particles.push_back(p);
        ccake::Particle<2>* pPtr = &sysPtr->particles[0];
        pPtr->thermo.s = 1.1;
        pPtr->thermo.rhoB = 1.2;
        pPtr->thermo.rhoQ = 1.3;
        pPtr->thermo.rhoS = 1.4;
    }

    return sysPtr;
}

void generate_fake_eos(){
    std::ofstream eos_file;
    ///Hypercube of points - Size 2x2x2x2 - 16 points
    double S_max = 2.0;
    double ds = 1;
    int Ns = 3;
    double rhoB_max = 2.0;
    double dB = 1;
    int NB = 3;
    double rhoQ_max = 2.0;
    double dQ = 1;
    int NQ = 3;
    double rhoS_max = 2.0;
    double dS = 1;
    int NS = 3;
    eos_file.open("eos_fake.dat", std::ios::out | std::ios::trunc);
    if (eos_file.is_open()){
      formatted_output::update("Writting header to 'eos_fake.dat'");
      eos_file << "#Smax= " << S_max << " ds= " << ds << " Ns= " << Ns << endl;
      eos_file << "#rhoBmax= " << rhoB_max << " drhoB= " << dB << " NB= " << NB << endl;
      eos_file << "#rhoQmax= " << rhoQ_max << " drhoQ= " << dQ << " NQ= " << NQ << endl;
      eos_file << "#rhoSmax= " << rhoS_max << " drhoS= " << dS << " NS= " << NS << endl;
      eos_file << "#s rhoB rhoQ rhoS e p cs2 muB muQ muS T s\n";
    } else {
      std::cerr << "File 'eos_fake.dat' could not be opened!\n";
      abort();
    }
    //Define our fake EoS
    auto fake_eos = [](double s, double B, double Q, double S){
        double e = s+B+Q+S;
        ///Return tuple with the thermodinamical quantities
        /// T muB muQ muS e p cs2 dwds dwdB dwdQ dwdS
        return std::make_tuple(e, e*2, e*3, e*4, e*5, e*6, e*7, e*8, e*9, e*10, e*11);
    };
    //Write eos to file
    for (int is=0; is<Ns; is++){
        for (int iB=0; iB<NB; iB++){
            for (int iQ=0; iQ<NQ; iQ++){
                for (int iS=0; iS<NS; iS++){
                    double s = is*ds;
                    double B = iB*dB;
                    double Q = iQ*dQ;
                    double S = iS*dS;
                    auto eos = fake_eos(s, B, Q, S);
                    eos_file << s << " " << B << " " << Q << " " << S << " ";
                    eos_file << std::get<0>(eos) << " " << std::get<1>(eos) << " " << std::get<2>(eos) << " " << std::get<3>(eos) << " ";
                    eos_file << std::get<4>(eos) << " " << std::get<5>(eos) << " " << std::get<6>(eos) << " " << std::get<7>(eos) << " ";
                    eos_file << std::get<8>(eos) << " " << std::get<9>(eos) << " " << std::get<10>(eos) << "fake_eos" << std::endl;
                }
            }
        }
    }
    eos_file.close();
}

TEST(eos_interpolator, linear_table){

    generate_fake_eos();
    std::shared_ptr<cc::SystemState<2>> systemPtr = prepare_system();
    cc::EoS_Interpolator eos("eos_fake.dat");

    double const expectedT = 5.;
    double const expectedmuB = 2*expectedT;
    double const expectedmuQ = 3*expectedT;
    double const expectedmuS = 4*expectedT;
    double const expectede = 5*expectedT;
    double const expectedp = 6*expectedT;
    double const expectedcs2 = 7*expectedT;
    double const expecteddwds = 8*expectedT;
    double const expecteddwdB = 9*expectedT;
    double const expecteddwdQ = 10*expectedT;
    double const expecteddwdS = 11*expectedT;


    systemPtr->allocate_cabana_particles();
    eos.fill_thermodynamics(systemPtr->cabana_particles,2);
    auto particles = Cabana::create_mirror_view_and_copy(Kokkos::HostSpace(), systemPtr->cabana_particles);
    CREATE_VIEW(host_,particles)

    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::T), expectedT);
    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::muB), expectedmuB);
    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::muQ), expectedmuQ);
    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::muS), expectedmuS);

    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::e), expectede);
    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::p), expectedp);
    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::cs2), expectedcs2);

    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::dwds), expecteddwds);
    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::dwdB), expecteddwdB);
    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::dwdQ), expecteddwdQ);
    ASSERT_DOUBLE_EQ(host_thermo(0,cc::thermo_info::dwdS), expecteddwdS);


}

int main(int argc, char* argv[]) {

    // Initialize Kokkos
    Kokkos::ScopeGuard guard(argc, argv);

    // Run the tests
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}