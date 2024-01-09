#include <memory>

#include <Cabana_Core.hpp>
#include "H5Cpp.h"
#include <omp.h>

#include "particle.h"
#include "eos.h"
#include "settings.h"
#include "system_state.h"
#include "thermodynamic_info.h"
#include "eos_interpolator.h"
#include "stopwatch.h"


//m#define ENABLE_TEST1
//#define ENABLE_TEST2
//#define SKIP_EOS_TABLE_GENERATION
#define SAVE_DOMAIN

int main(int argc, char** argv){

    //Get number of threads
    const int nthreads = omp_get_max_threads();

    ccake::EquationOfState eos[nthreads];
    int idx_a[nthreads];
    std::shared_ptr<Settings> settings = std::make_shared<Settings>();

    //Setup EoS
    //settings->eos_type = "conformal"; // "table", "tanh_conformal", "conformal" or "conformal_diagonal"
    settings->eos_type = "table"; // "table", "tanh_conformal", "conformal" or "conformal_diagonal"
    settings->eos_path = "EoS/Houston"; //Only if table

    //Cannot be initialized in parallel to avoid simultaneous access to HDF5 file
    for(int i=0; i<nthreads; ++i){
      eos[i].set_SettingsPtr(settings);
      eos[i].init();
      idx_a[i] = 0;
    }


    //Setup grid
    const double s_min = 0.0001;
    const double B_min = 0.0001;
    const double S_min = 0.0001;
    const double Q_min = 0.0001;

    //Limits obtained by looking at a sample of 100 iccing events
    const double s_max = 70.0;
    const double B_max = 5.0;
    const double S_max = 7.0;
    const double Q_max = 8.6;
    //Values suitable for conformal EoS
    //const int Ns = 40000+1;
    //const int NB = 1+1;
    //const int NQ = 1+1;
    //const int NS = 1+1;
    const int Ns = 31+1;
    const int NB = 31+1;
    const int NQ = 31+1;
    const int NS = 31+1;
    const double ds = (s_max-s_min)/(Ns-1);
    const double dB = (B_max-B_min)/(NB-1);
    const double dS = (S_max-s_min)/(NS-1);
    const double dQ = (Q_max-Q_min)/(NQ-1);

    //Estimate final size of EoS table
    const long Nsize = Ns*NB*NS*NQ;
    const long Nbytes = Nsize*sizeof(double)*11; //11 thermodynamic quantities stored
    const double NbytesGB = Nbytes/1024./1024./1024.;
    cout << "EoS table size: " << std::setprecision(4) << NbytesGB << " GB" << endl;
    int idx=0;
    #ifndef SKIP_EOS_TABLE_GENERATION
    //Write header to HDF5 file
    H5::H5File eos_h5_file("eos_table.h5", H5F_ACC_TRUNC); //Open file in binary mode
    if (eos_h5_file.getId() < 0){
      std::cerr << "File 'eos_table.h5' could not be opened!\n";
      abort();
    }

    //Create datasets
    hsize_t dims[4] = {Ns, NB, NQ, NS};
    H5::DataSpace table_dataspace(4, dims);
    #ifdef SAVE_DOMAIN
    H5::DataSet s_dataset = eos_h5_file.createDataSet("s", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet rhoB_dataset = eos_h5_file.createDataSet("rhoB", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet rhoQ_dataset = eos_h5_file.createDataSet("rhoQ", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet rhoS_dataset = eos_h5_file.createDataSet("rhoS", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    #endif
    H5::DataSet T_dataset = eos_h5_file.createDataSet("T", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet muB_dataset = eos_h5_file.createDataSet("muB", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet muQ_dataset = eos_h5_file.createDataSet("muQ", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet muS_dataset = eos_h5_file.createDataSet("muS", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet e_dataset = eos_h5_file.createDataSet("e", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet p_dataset = eos_h5_file.createDataSet("p", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet cs2_dataset = eos_h5_file.createDataSet("cs2", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet dwds_dataset = eos_h5_file.createDataSet("dwds", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet dwdB_dataset = eos_h5_file.createDataSet("dwdB", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet dwdQ_dataset = eos_h5_file.createDataSet("dwdQ", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet dwdS_dataset = eos_h5_file.createDataSet("dwdS", H5::PredType::NATIVE_DOUBLE, table_dataspace);
    H5::DataSet eos_type_dataset = eos_h5_file.createDataSet("eos_type", H5::PredType::NATIVE_INT, table_dataspace);

    ///Create attributes
    H5::DataSpace scalar_type(H5S_SCALAR);
    H5::Attribute smin_attr = eos_h5_file.createAttribute("s_min", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute smax_attr = eos_h5_file.createAttribute("s_max", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute ds_attr = eos_h5_file.createAttribute("ds", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute rhoBmin_attr = eos_h5_file.createAttribute("rhoB_min", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute rhoBmax_attr = eos_h5_file.createAttribute("rhoB_max", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute drhoB_attr = eos_h5_file.createAttribute("drhoB", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute rhoSmin_attr = eos_h5_file.createAttribute("rhoS_min", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute rhoSmax_attr = eos_h5_file.createAttribute("rhoS_max", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute drhoS_attr = eos_h5_file.createAttribute("drhoS", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute rhoQmin_attr = eos_h5_file.createAttribute("rhoQ_min", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute rhoQmax_attr = eos_h5_file.createAttribute("rhoQ_max", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute drhoQ_attr = eos_h5_file.createAttribute("drhoQ", H5::PredType::NATIVE_DOUBLE, scalar_type);
    H5::Attribute sAxis_attr = eos_h5_file.createAttribute("s_Axis", H5::PredType::NATIVE_INT, scalar_type);
    H5::Attribute rhoBAxis_attr = eos_h5_file.createAttribute("rhoB_Axis", H5::PredType::NATIVE_INT, scalar_type);
    H5::Attribute rhoSAxis_attr = eos_h5_file.createAttribute("rhoS_Axis", H5::PredType::NATIVE_INT, scalar_type);
    H5::Attribute rhoQAxis_attr = eos_h5_file.createAttribute("rhoQ_Axis", H5::PredType::NATIVE_INT, scalar_type);

    hsize_t dims_eos_types[1] = {4}; //4 known EoS types
    H5::DataSpace eos_types_dataspace(1, dims_eos_types);
    H5::StrType strtype(H5::PredType::C_S1, H5T_VARIABLE);
    H5::Attribute eos_types_codes = eos_h5_file.createAttribute("eos_type_codes", H5::PredType::NATIVE_INT, eos_types_dataspace);
    H5::Attribute eos_names_attr = eos_h5_file.createAttribute("eos_type_names", strtype, eos_types_dataspace);

    //Write attributes
    smin_attr.write(H5::PredType::NATIVE_DOUBLE, &s_min);
    smax_attr.write(H5::PredType::NATIVE_DOUBLE, &s_max);
    ds_attr.write(H5::PredType::NATIVE_DOUBLE, &ds);
    rhoBmin_attr.write(H5::PredType::NATIVE_DOUBLE, &B_min);
    rhoBmax_attr.write(H5::PredType::NATIVE_DOUBLE, &B_max);
    drhoB_attr.write(H5::PredType::NATIVE_DOUBLE, &dB);
    rhoSmin_attr.write(H5::PredType::NATIVE_DOUBLE, &S_min);
    rhoSmax_attr.write(H5::PredType::NATIVE_DOUBLE, &S_max);
    drhoS_attr.write(H5::PredType::NATIVE_DOUBLE, &dS);
    rhoQmin_attr.write(H5::PredType::NATIVE_DOUBLE, &Q_min);
    rhoQmax_attr.write(H5::PredType::NATIVE_DOUBLE, &Q_max);
    drhoQ_attr.write(H5::PredType::NATIVE_DOUBLE, &dQ);

    int axis=0;
    sAxis_attr.write(H5::PredType::NATIVE_INT, &axis);
    axis=1;
    rhoBAxis_attr.write(H5::PredType::NATIVE_INT, &axis);
    axis=2;
    rhoSAxis_attr.write(H5::PredType::NATIVE_INT, &axis);
    axis=3;
    rhoQAxis_attr.write(H5::PredType::NATIVE_INT, &axis);

    //int eos_codes[4] = {0, 1, 2, 3};
    //eos_names_attr.write(H5::PredType::NATIVE_INT, eos_codes);
    //std::string eos_names[4] = {"table", "tanh_conformal", "conformal", "conformal_diagonal"};
    //eos_names_attr.write(strtype, eos_names);

    //Storage for EoS table
    const int N = Ns*NB*NS*NQ;
    std::array<double, N>* s_array = new std::array<double, N>; //Let us manage this, since we want to free it as soon as possible
    std::array<double, N>* rhoB_array = new std::array<double, N>;
    std::array<double, N>* rhoQ_array = new std::array<double, N>;
    std::array<double, N>* rhoS_array = new std::array<double, N>;
    std::array<double, N>* T_array = new std::array<double, N>;
    std::array<double, N>* muB_array = new std::array<double, N>;
    std::array<double, N>* muQ_array = new std::array<double, N>;
    std::array<double, N>* muS_array = new std::array<double, N>;
    std::array<double, N>* e_array = new std::array<double, N>;
    std::array<double, N>* p_array = new std::array<double, N>;
    std::array<double, N>* cs2_array = new std::array<double, N>;
    std::array<double, N>* dwds_array = new std::array<double, N>;
    std::array<double, N>* dwdB_array = new std::array<double, N>;
    std::array<double, N>* dwdS_array = new std::array<double, N>;
    std::array<double, N>* dwdQ_array = new std::array<double, N>;

    //Get start time
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    int ith=0;
    #pragma omp parallel for collapse(4) schedule(dynamic) private(ith,idx)
    for (int is=0; is<Ns;++is)
    for (int iB=0; iB<NB;++iB)
    for (int iS=0; iS<NS;++iS)
    for (int iQ=0; iQ<NQ;++iQ)
    {

        ith = omp_get_thread_num();
        idx = is*NB*NQ*NS + iB*NQ*NS + iS*NQ + iQ;

        double s_In = s_min+is*ds;
        double rhoB_In = B_min+iB*dB;
        double rhoS_In = S_min+iS*dS;
        double rhoQ_In = Q_min+iQ*dQ;
        double T = 750./197.;
        double muB = 0.0;
        double muQ = 0.0;
        double muS = 0.0;

        #ifdef SAVE_DOMAIN
        (*s_array)[idx] = s_In;
        (*rhoB_array)[idx] = rhoB_In;
        (*rhoQ_array)[idx] = rhoQ_In;
        (*rhoS_array)[idx] = rhoS_In;
        #endif
        eos[ith].tbqs(T, muB, muQ, muS, "default");
        bool update_s_success = eos[ith].update_s(s_In, rhoB_In, rhoS_In, rhoQ_In, false);
        ccake::thermodynamic_info thermo;
        if ( update_s_success ){
          eos[ith].set_thermo( thermo );
        } else {
            std::cerr << "EoS root-finding failed at grid point: " << is << " " << iB << " " << iQ << " " << iS << std::endl;
            abort();
        }
        (*T_array)[idx] = thermo.T;
        (*muB_array)[idx] = thermo.muB;
        (*muQ_array)[idx] = thermo.muQ;
        (*muS_array)[idx] = thermo.muS;
        (*e_array)[idx] = thermo.e;
        (*p_array)[idx] = thermo.p;
        (*cs2_array)[idx] = thermo.cs2;
        (*dwds_array)[idx] = thermo.dwds;
        (*dwdB_array)[idx] = thermo.dwdB;
        (*dwdS_array)[idx] = thermo.dwdS;
        (*dwdQ_array)[idx] = thermo.dwdQ;

        //Get end time
        if(ith==0 && idx_a[0]%100==0){
            #pragma omp critical
            {
              int idx_s=0;
              for (int i=1; i<nthreads; ++i) idx_s += idx_a[i];
              end = std::chrono::system_clock::now();
              std::chrono::duration<double> elapsed_seconds = end-start;
              std::cout << "EoS table generation: " << idx_s << " / " << N << " = " << (double)idx_s/N*100. << "% done. Elapsed time: " << elapsed_seconds.count() << "s";
              std::cout << " Estimated time remaining: " << elapsed_seconds.count()/(double)idx_s*(N-idx_s) << "s" << std::endl;
            }
        }
        idx_a[ith]++;
    }
    end = std::chrono::system_clock::now();
    cout << "EoS table generation: 100% done. Time elapsed: " << std::chrono::duration<double>(end-start).count() << "s" << endl;


    //Write to file
    #ifdef SAVE_DOMAIN
    s_dataset.write(s_array->data(), H5::PredType::NATIVE_DOUBLE);
    rhoB_dataset.write(rhoB_array->data(), H5::PredType::NATIVE_DOUBLE);
    rhoQ_dataset.write(rhoQ_array->data(), H5::PredType::NATIVE_DOUBLE);
    rhoS_dataset.write(rhoS_array->data(), H5::PredType::NATIVE_DOUBLE);
    #endif
    T_dataset.write(T_array->data(), H5::PredType::NATIVE_DOUBLE);
    muB_dataset.write(muB_array->data(), H5::PredType::NATIVE_DOUBLE);
    muQ_dataset.write(muQ_array->data(), H5::PredType::NATIVE_DOUBLE);
    muS_dataset.write(muS_array->data(), H5::PredType::NATIVE_DOUBLE);
    e_dataset.write(e_array->data(), H5::PredType::NATIVE_DOUBLE);
    p_dataset.write(p_array->data(), H5::PredType::NATIVE_DOUBLE);
    cs2_dataset.write(cs2_array->data(), H5::PredType::NATIVE_DOUBLE);
    dwds_dataset.write(dwds_array->data(), H5::PredType::NATIVE_DOUBLE);
    dwdB_dataset.write(dwdB_array->data(), H5::PredType::NATIVE_DOUBLE);
    dwdS_dataset.write(dwdS_array->data(), H5::PredType::NATIVE_DOUBLE);
    dwdQ_dataset.write(dwdQ_array->data(), H5::PredType::NATIVE_DOUBLE);

    delete s_array, rhoB_array, rhoQ_array, rhoS_array,
           T_array, muB_array, muQ_array, muS_array,
           e_array, p_array, cs2_array,
           dwds_array, dwdB_array, dwdQ_array, dwdS_array; //Free memory
    eos_h5_file.close();

    #endif

    //Perform a few tests to ensure that the EoS will be interpolated correctly
    Kokkos::ScopeGuard scope_guard(argc, argv); //Create Kokkos scope guard to
    // ensure that Kokkos is properly initialized and finalized
    ccake::EoS_Interpolator eos_interpolator("eos_table.h5");
    #ifdef ENABLE_TEST1
    std::shared_ptr<ccake::SystemState<1>> systemPtr = std::make_shared<ccake::SystemState<1>>(settings);
    //Test 1: Check that the EoS is correctly interpolated at the grid points
    systemPtr->n_particles = Ns*NB*NQ*NS;
    systemPtr->t = 2.;
    for (int is=0; is<Ns;++is)
    for (int iB=0; iB<NB;++iB)
    for (int iQ=0; iQ<NQ;++iQ)
    for (int iS=0; iS<NS;++iS){
      ccake::Particle<1> particle;
      systemPtr->particles.push_back(particle);
    }


    settings->xmin = -4;
    settings->ymin = -6;
    settings->etamin = -2;
    settings->hT = 1;
    settings->using_shear = true;
    settings->t0 = 2;
    systemPtr->initialize();


    idx = 0;
    for (int is=0; is<Ns;++is)
    for (int iB=0; iB<NB;++iB)
    for (int iQ=0; iQ<NQ;++iQ)
    for (int iS=0; iS<NS;++iS){
      ccake::Particle<1>* p = &systemPtr->particles[idx];
      p->ID = idx;
      p->thermo.s = is*ds;
      p->thermo.rhoB = iB*dB;
      p->thermo.rhoQ = iQ*dQ;
      p->thermo.rhoS = iS*dS;
      idx++;
    }
    systemPtr->allocate_cabana_particles();
    eos_interpolator.fill_thermodynamics(systemPtr->cabana_particles, 0.0);
    systemPtr->copy_device_to_host();


    idx = 0;
    cout << "Testing EoS interpolation at grid points" << endl;
    for (int is=0; is<Ns;++is)
    for (int iB=0; iB<NB;++iB)
    for (int iQ=0; iQ<NQ;++iQ)
    for (int iS=0; iS<NS;++iS){

      eos[0].tbqs(800./197., .0, .0, .0, "default");
      double s_In = is*ds;
      double rhoB_In = iB*dB;
      double rhoQ_In = iQ*dQ;
      double rhoS_In = iS*dS;
      bool update_s_success = eos[0].update_s(s_In, rhoB_In, rhoS_In, rhoQ_In, false);

      ccake::thermodynamic_info thermo;
      if ( update_s_success ) eos[0].set_thermo( thermo );
      //Check that the EoS is correctly interpolated
      if ( std::abs(thermo.T - systemPtr->particles[idx].thermo.T ) > 1.e-13 ){
       std::cerr << "EoS interpolation failed at point = " << idx << std::endl;
       std::cerr << "s_rootfinder      = " << std::setprecision(13) << thermo.s << " "
                 << "rhoB_rootfinder   = " << std::setprecision(13) << thermo.rhoB << " "
                 << "rhoQ_rootfinder   = " << std::setprecision(13) << thermo.rhoQ << " "
                 << "rhoS_rootfinder   = " << std::setprecision(13) << thermo.rhoS << " "
                 << "T_rootfinder      = " << std::setprecision(13) << thermo.T << " " << endl;
        std::cerr << "s_interpolator    = " << std::setprecision(13) << systemPtr->particles[idx].thermo.s << " "
                  << "rhoB_interpolator = " << std::setprecision(13) << systemPtr->particles[idx].thermo.rhoB << " "
                  << "rhoQ_interpolator = " << std::setprecision(13) << systemPtr->particles[idx].thermo.rhoQ << " "
                  << "rhoS_interpolator = " << std::setprecision(13) << systemPtr->particles[idx].thermo.rhoS << " "
                  << "T_interpolator    = " << std::setprecision(13) << systemPtr->particles[idx].thermo.T << " " << endl;

       abort();
      }

      idx++;
    }
    cout << "Test 1 passed" << endl;
    #endif

    #ifdef ENABLE_TEST2
    //Test 2: Check that the EoS is correctly interpolated at random points
    cout << "Testing EoS interpolation at random points" << endl;

    std::shared_ptr<ccake::SystemState<1>> systemPtr2 = std::make_shared<ccake::SystemState<1>>(settings);
    //Test 1: Check that the EoS is correctly interpolated at the grid points
    systemPtr2->n_particles = 1000;
    systemPtr2->t = 2.;
    for (int i=0; i<1000;++i)
    {
      ccake::Particle<1> particle;
      systemPtr2->particles.push_back(particle);
    }

    settings->xmin = -4;
    settings->ymin = -6;
    settings->etamin = -2;
    settings->hT = 1;
    settings->using_shear = true;
    settings->t0 = 2;
    systemPtr2->initialize();


    //set random seed
    srand(42); //use a fixed seed for reproducibility
               //42, since it is the answer to the ultimate question of life, the universe and everything
               //and so, it is also the answer for checking the EoS interpolation ;)

    double s_In[1000], rhoB_In[1000], rhoQ_In[1000], rhoS_In[1000];
    for (int i=0; i<1000; ++i){
      s_In[i] = (double)rand()/RAND_MAX*s_max;
      rhoB_In[i] = (double)rand()/RAND_MAX*B_max;
      rhoQ_In[i] = (double)rand()/RAND_MAX*Q_max;
      rhoS_In[i] = (double)rand()/RAND_MAX*S_max;
      ccake::Particle<1>* p = &systemPtr2->particles[i];
      p->ID = i;
      p->thermo.s = s_In[i];
      p->thermo.rhoB = rhoB_In[i];
      p->thermo.rhoQ = rhoQ_In[i];
      p->thermo.rhoS = rhoS_In[i];
    }
    systemPtr2->allocate_cabana_particles();
    eos_interpolator.fill_thermodynamics(systemPtr2->cabana_particles, 0.0);
    systemPtr2->copy_device_to_host();

    //Reset random seed
    srand(42);
    for (int i=0; i<1000; ++i){
      double const tol = 2.e-4;
      eos[0].tbqs(800./197., .0, .0, .0, "default");
      bool update_s_success = eos[0].update_s(s_In[i], rhoB_In[i], rhoS_In[i], rhoQ_In[i], false);

      ccake::thermodynamic_info thermo;
      if ( update_s_success ) eos[0].set_thermo( thermo );

      //Check that the EoS is correctly interpolated
      double deltaT = std::abs(thermo.T - systemPtr2->particles[i].thermo.T )/thermo.T;
      double deltaMuB = std::abs(thermo.muB - systemPtr2->particles[i].thermo.muB )/thermo.muB;
      double deltaMuQ = std::abs(thermo.muQ - systemPtr2->particles[i].thermo.muQ )/thermo.muQ;
      double deltaMuS = std::abs(thermo.muS - systemPtr2->particles[i].thermo.muS )/thermo.muS;

      double deltaE = std::abs(thermo.e - systemPtr2->particles[i].thermo.e )/thermo.e;
      double deltaP = std::abs(thermo.p - systemPtr2->particles[i].thermo.p )/thermo.p;
      double deltaCs2 = std::abs(thermo.cs2 - systemPtr2->particles[i].thermo.cs2 )/thermo.cs2;

      double deltaDwds = std::abs(thermo.dwds - systemPtr2->particles[i].thermo.dwds )/thermo.dwds;
      double deltaDwdB = std::abs(thermo.dwdB - systemPtr2->particles[i].thermo.dwdB )/thermo.dwdB;
      double deltaDwdQ = std::abs(thermo.dwdQ - systemPtr2->particles[i].thermo.dwdQ )/thermo.dwdQ;
      double deltaDwdS = std::abs(thermo.dwdS - systemPtr2->particles[i].thermo.dwdS )/thermo.dwdS;

      if (deltaT > tol)
        std::cerr << "T interpolation failed at point = (" << s_In[i] << ", " << rhoB_In[i] << ", " << rhoQ_In[i] << ", " << rhoS_In[i] << ") "
          << " Delta = "  << deltaT << " T = " << thermo.T << std::endl;
      if (deltaMuB > tol)
        std::cerr << "muB interpolation failed at point = (" << s_In[i] << ", " << rhoB_In[i] << ", " << rhoQ_In[i] << ", " << rhoS_In[i] << ") "
          << " Delta = "  << deltaMuB
                  << " muB = " << thermo.muB << " muB_interp = " << systemPtr2->particles[i].thermo.muB << std::endl;
      if (deltaMuQ > tol)
        std::cerr << "muQ interpolation failed at point = (" << s_In[i] << ", " << rhoB_In[i] << ", " << rhoQ_In[i] << ", " << rhoS_In[i] << ") "
          << " Delta = "  << deltaMuQ
                  << " muQ = " << thermo.muQ << " muQ_interp = " << systemPtr2->particles[i].thermo.muQ << std::endl;
      if (deltaMuS > tol)
        std::cerr << "muS interpolation failed at point = (" << s_In[i] << ", " << rhoB_In[i] << ", " << rhoQ_In[i] << ", " << rhoS_In[i] << ") "
          << " Delta = "  << deltaMuS
                  << " muS = " << thermo.muS << " muS_interp = " << systemPtr2->particles[i].thermo.muS << std::endl;
      if (deltaE > tol)
        std::cerr << "e interpolation failed at point = " << i << " Delta = "  << deltaE << std::endl;
      if (deltaP > tol)
        std::cerr << "p interpolation failed at point = " << i << " Delta = "  << deltaP << std::endl;
      if (deltaCs2 > tol)
        std::cerr << "cs2 interpolation failed at point = " << i << " Delta = "  << deltaCs2 << std::endl;
      if (deltaDwds > tol)
        std::cerr << "dwds interpolation failed at point = " << i << " Delta = "  << deltaDwds << std::endl;
      if (deltaDwdB > tol)
        std::cerr << "dwdB interpolation failed at point = " << i << " Delta = "  << deltaDwdB << std::endl;
      if (deltaDwdQ > tol)
        std::cerr << "dwdQ interpolation failed at point = " << i << " Delta = "  << deltaDwdQ << std::endl;
      if (deltaDwdS > tol)
        std::cerr << "dwdS interpolation failed at point = " << i << " Delta = "  << deltaDwdS << std::endl;

    }
    cout << "Test 2 finished." << endl;
    #endif
}