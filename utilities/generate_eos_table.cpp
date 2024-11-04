/// @file generate_eos_table
/// @author Willian M. Serenone
/// @date 2024-03-19
/// @brief Generate a table of EoS values for a given EoS
/// @details This program generates a table of EoS values for a given EoS. The
/// table is saved in a HDF5 file, and the table is divided in subtables, each
/// corresponding to a decade in each of the 4 dimensions. The table is
/// generated in parallel, and the number of threads is set by the user.


#include <memory>
#include <filesystem>

#include "H5Cpp.h"
#include <omp.h>

#include "eos.h"
#include "settings.h"
#include "thermodynamic_info.h"
#include "stopwatch.h"
#include <assert.h>

#define SAVE_DOMAIN

namespace fs = std::filesystem;

int check_args(int argc, char** argv){
  //Usage: ./generate_eos_table <nthreads> <eos_type> <eos_path> <output_file>
  if (argc != 5){
    std::cerr << "Usage: " << argv[0] <<
    " <nthreads> <eos_type> <eos_path> <output_file>" << std::endl;
    std::cerr << "nthreads: number of threads to use" << std::endl;
    std::cerr << "eos_type: 'table', 'tanh_conformal', 'conformal' or 'conformal_diagonal'" << std::endl;
    std::cerr << "eos_path: path to EoS table" << std::endl;
    std::cerr << "output_file: name of the output file" << std::endl;
    return 1;
  } else {
    std::cout << "Number of threads: " << argv[1] << std::endl;
    std::cout << "EoS type.........: " << argv[2] << std::endl;
    std::cout << "EoS path.........: " << argv[3] << std::endl;
    std::cout << "Output file......: " << argv[4] << std::endl;
    return 0;
  }
}

class Inverter{
  private:
    double x_min[4]; //max values for s, rhoB, rhoQ, rhoS
    double x_max[4]; //min values for s, rhoB, rhoQ, rhoS
    int N[4]; //Number of points in each dimension
    double dx[4]; //step size in each dimension

    //HDF5 file
    fs::path output_file;
    H5::H5File eos_h5_file;

    ///Possible datasets
    #ifdef SAVE_DOMAIN
    H5::DataSet s_dataset;
    H5::DataSet rhoB_dataset;
    H5::DataSet rhoQ_dataset;
    H5::DataSet rhoS_dataset;
    #endif
    H5::DataSet T_dataset;
    H5::DataSet muB_dataset;
    H5::DataSet muQ_dataset;
    H5::DataSet muS_dataset;
    H5::DataSet e_dataset;
    H5::DataSet p_dataset;
    H5::DataSet cs2_dataset;
    H5::DataSet dwds_dataset;
    H5::DataSet dwdB_dataset;
    H5::DataSet dwdQ_dataset;
    H5::DataSet dwdS_dataset;
    H5::DataSet eos_type_dataset;

    //Auxiliary datamembers
    ccake::EquationOfState* eos;
    std::shared_ptr<Settings> settings;
    int nthreads;

    bool eos_initialized = false;
    bool range_set = false;

    ///@brief Helper function to create datasets in the HDF5 group
    void create_datasets(H5::Group group){
      hsize_t dims[4] = {N[0], N[1], N[2], N[3]};
      H5::DataSpace table_dataspace(4, dims);
      #ifdef SAVE_DOMAIN
      s_dataset = group.createDataSet("s", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      rhoB_dataset = group.createDataSet("rhoB", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      rhoQ_dataset = group.createDataSet("rhoQ", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      rhoS_dataset = group.createDataSet("rhoS", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      #endif
      T_dataset = group.createDataSet("T", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      muB_dataset = group.createDataSet("muB", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      muQ_dataset = group.createDataSet("muQ", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      muS_dataset = group.createDataSet("muS", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      e_dataset = group.createDataSet("e", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      p_dataset = group.createDataSet("p", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      cs2_dataset = group.createDataSet("cs2", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      dwds_dataset = group.createDataSet("dwds", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      dwdB_dataset = group.createDataSet("dwdB", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      dwdQ_dataset = group.createDataSet("dwdQ", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      dwdS_dataset = group.createDataSet("dwdS", H5::PredType::NATIVE_DOUBLE, table_dataspace);
      eos_type_dataset = group.createDataSet("eos_type", H5::PredType::NATIVE_INT, table_dataspace);
    }

    void create_header(H5::Group group){
      H5::DataSpace scalar_type(H5S_SCALAR);
      H5::Attribute smin_attr = group.createAttribute("s_min", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute smax_attr = group.createAttribute("s_max", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute ds_attr = group.createAttribute("ds", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute rhoBmin_attr = group.createAttribute("rhoB_min", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute rhoBmax_attr = group.createAttribute("rhoB_max", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute drhoB_attr = group.createAttribute("drhoB", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute rhoSmin_attr = group.createAttribute("rhoS_min", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute rhoSmax_attr = group.createAttribute("rhoS_max", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute drhoS_attr = group.createAttribute("drhoS", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute rhoQmin_attr = group.createAttribute("rhoQ_min", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute rhoQmax_attr = group.createAttribute("rhoQ_max", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute drhoQ_attr = group.createAttribute("drhoQ", H5::PredType::NATIVE_DOUBLE, scalar_type);
      H5::Attribute sAxis_attr = group.createAttribute("s_Axis", H5::PredType::NATIVE_INT, scalar_type);
      H5::Attribute rhoBAxis_attr = group.createAttribute("rhoB_Axis", H5::PredType::NATIVE_INT, scalar_type);
      H5::Attribute rhoSAxis_attr = group.createAttribute("rhoS_Axis", H5::PredType::NATIVE_INT, scalar_type);
      H5::Attribute rhoQAxis_attr = group.createAttribute("rhoQ_Axis", H5::PredType::NATIVE_INT, scalar_type);

      smin_attr.write(H5::PredType::NATIVE_DOUBLE, &x_min[0]);
      smax_attr.write(H5::PredType::NATIVE_DOUBLE, &x_max[0]);
      ds_attr.write(H5::PredType::NATIVE_DOUBLE, &dx[0]);
      rhoBmin_attr.write(H5::PredType::NATIVE_DOUBLE, &x_min[1]);
      rhoBmax_attr.write(H5::PredType::NATIVE_DOUBLE, &x_max[1]);
      drhoB_attr.write(H5::PredType::NATIVE_DOUBLE, &dx[1]);
      rhoSmin_attr.write(H5::PredType::NATIVE_DOUBLE, &x_min[2]);
      rhoSmax_attr.write(H5::PredType::NATIVE_DOUBLE, &x_max[2]);
      drhoS_attr.write(H5::PredType::NATIVE_DOUBLE, &dx[2]);
      rhoQmin_attr.write(H5::PredType::NATIVE_DOUBLE, &x_min[3]);
      rhoQmax_attr.write(H5::PredType::NATIVE_DOUBLE, &x_max[3]);
      drhoQ_attr.write(H5::PredType::NATIVE_DOUBLE, &dx[3]);

      int axis=0;
      sAxis_attr.write(H5::PredType::NATIVE_INT, &axis);
      axis=1;
      rhoBAxis_attr.write(H5::PredType::NATIVE_INT, &axis);
      axis=2;
      rhoSAxis_attr.write(H5::PredType::NATIVE_INT, &axis);
      axis=3;
      rhoQAxis_attr.write(H5::PredType::NATIVE_INT, &axis);
    }

  public:
    Inverter() = delete;
    /// @brief Create an inverter object
    /// @details The creator will be responsible for creating the HDF5. It will
    /// not overwrite an existing file.
    /// @param output The path where the inverted table will be saved
    Inverter(fs::path output, int nthreads_in) :
      output_file(output),  nthreads(nthreads_in)
    {

      output_file = output; //Path to output file
      //Check if the file exists and abort if it does
      if (fs::exists(output_file)){
        std::cerr << "Output file " << output_file << " already exists!" << std::endl;
        abort();
      }

      eos_h5_file = H5::H5File(output_file.string() , H5F_ACC_TRUNC);
      //Abort if the file could not be opened
      if (eos_h5_file.getId() < 0){
        std::cerr << "File 'eos_table.h5' could not be opened!\n";
        abort();
      }
    }

    void init_ccake(std::string eos_type, std::string eos_path){
      settings = std::make_shared<Settings>();
      settings->eos_type = eos_type;
      settings->eos_path = eos_path;
      eos = new ccake::EquationOfState[nthreads];
      for(int i=0; i<nthreads; ++i){
        eos[i].set_SettingsPtr(settings);
        eos[i].init();
      }
      eos_initialized = true;
    }

    /// @brief Set the range of the table
    /// @details The table always cover a full decade in each axis. We specify
    /// which decade as an integer, which is the exponent of the minimum value.
    /// E.g., a value of -2 means that the minimum value is 0.01, and the
    /// maximum value is 0.1. The number of points in each dimension is always
    /// 10. The step size is calculated as max/10.
    /// @param exponent An array of 4 integers, each representing the exponent
    /// of the range in one of the 4 directions.
    void set_range(int exponent[4]){
      for (int j=0; j<4; ++j){
        x_min[j] = pow(10., exponent[j]);
        x_max[j] = pow(10., exponent[j]+1);
        N[j] = 20;
        dx[j] = x_max[j]/N[j];
      }
      range_set = true;
      std::cout << "Range for s set to...: " << x_min[0] << " " << x_max[0] << " " << N[0] << " " << dx[0] << std::endl;
      std::cout << "Range for rhoB set to: " << x_min[1] << " " << x_max[1] << " " << N[1] << " " << dx[1] << std::endl;
      std::cout << "Range for rhoS set to: " << x_min[2] << " " << x_max[2] << " " << N[2] << " " << dx[2] << std::endl;
      std::cout << "Range for rhoQ set to: " << x_min[3] << " " << x_max[3] << " " << N[3] << " " << dx[3] << std::endl;
    }

    void set_range(std::initializer_list<int> exponent_in){
      assert(exponent_in.size() == 4);
      int exponent[4];
      std::copy(exponent_in.begin(), exponent_in.end(), exponent);
      set_range(exponent);
    }

    //Getters
    void get_xmin(int *x){for (int j=0; j<4; ++j) x[j] = x_min[j];}
    void get_xmax(int *x){for (int j=0; j<4; ++j) x[j] = x_max[j];}
    void get_N(int *x){for (int j=0; j<4; ++j) x[j] = N[j];}
    void get_dx(double *x){for (int j=0; j<4; ++j) x[j] = dx[j];}

    /// @brief Find the table, respecting the limits set in the class
    void create_table(std::string folder_name){

      assert(range_set);
      assert(eos_initialized);

      //Get table limits
      int xmin[4]; get_xmin(xmin);
      int N[4]; get_N(N);
      double dx[4]; get_dx(dx);
      int Nsize = N[0]*N[1]*N[2]*N[3];

      //Create a folder inside the hdf5 file
      cout << "Creating table " << "/"+folder_name << endl;
      H5::Group group = eos_h5_file.createGroup("/"+folder_name);
      create_datasets(group);
      create_header(group);

    //==========================================================================
    //Create memory storage for results
    //==========================================================================
    std::size_t tot_size = N[0]*N[1]*N[2]*N[3];
    //Create at heap to avoid stack overflow and to release as soon as possible
    #ifdef SAVE_DOMAIN
    std::vector<double>* s_array = new std::vector<double>(tot_size);
    std::vector<double>* rhoB_array = new std::vector<double>(tot_size);
    std::vector<double>* rhoQ_array = new std::vector<double>(tot_size);
    std::vector<double>* rhoS_array = new std::vector<double>(tot_size);
    #endif
    std::vector<double>* T_array = new std::vector<double>(tot_size);
    std::vector<double>* muB_array = new std::vector<double>(tot_size);
    std::vector<double>* muQ_array = new std::vector<double>(tot_size);
    std::vector<double>* muS_array = new std::vector<double>(tot_size);
    std::vector<double>* e_array = new std::vector<double>(tot_size);
    std::vector<double>* p_array = new std::vector<double>(tot_size);
    std::vector<double>* cs2_array = new std::vector<double>(tot_size);
    std::vector<double>* dwds_array = new std::vector<double>(tot_size);
    std::vector<double>* dwdB_array = new std::vector<double>(tot_size);
    std::vector<double>* dwdS_array = new std::vector<double>(tot_size);
    std::vector<double>* dwdQ_array = new std::vector<double>(tot_size);

    //==========================================================================
    //Main loop
    //==========================================================================

    //Get start time
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    int Ns = N[0];
    int NB = N[1];
    int NS = N[2];
    int NQ = N[3];
    double s_min = x_min[0];
    double B_min = x_min[1];
    double S_min = x_min[2];
    double Q_min = x_min[3];
    double ds = dx[0];
    double dB = dx[1];
    double dS = dx[2];
    double dQ = dx[3];

    std::cout << "EoS table generation parameters: "<< std::endl;
    std::cout << "Size..........: " << Ns << " " << NB << " "
                                    << NS << " " << NQ << std::endl;
    std::cout << "Minimum values: " << s_min << " " << B_min << " "
                                    << S_min << " " << Q_min << " " << std::endl;
    std::cout << "Step size.....: " << ds << " " << dB << " "
                                    << dS << " " << dQ << endl;

    const int Ntotal = Ns*NB*NS*NQ;
    int idx_a[nthreads]; for(int i=0; i<nthreads; ++i) idx_a[i]=0;
    #pragma omp parallel for collapse(4) schedule(dynamic)
    for (int is=0; is<Ns;++is)
    for (int iB=0; iB<NB;++iB)
    for (int iS=0; iS<NS;++iS)
    for (int iQ=0; iQ<NQ;++iQ)
    {

        int ith = omp_get_thread_num();
        int idx = is*NB*NQ*NS + iB*NQ*NS + iS*NQ + iQ;

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
        (*rhoS_array)[idx] = rhoS_In;
        (*rhoQ_array)[idx] = rhoQ_In;
        #endif
        eos[ith].tbqs(T, muB, muQ, muS, "default");
        bool update_s_success = eos[ith].update_s(s_In, rhoB_In,
                                                  rhoS_In, rhoQ_In, false);
        ccake::thermodynamic_info thermo;
        if ( update_s_success ){
          eos[ith].set_thermo( thermo );
        } else {
            std::cerr << "EoS root-finding failed at grid point: " << is << " " << iB << " " << iQ << " " << iS << std::endl;
            abort();
        }
        (*T_array)[idx] = thermo.T;
        (*muB_array)[idx] = thermo.muB;
        (*muS_array)[idx] = thermo.muS;
        (*muQ_array)[idx] = thermo.muQ;
        (*e_array)[idx] = thermo.e;
        (*p_array)[idx] = thermo.p;
        (*cs2_array)[idx] = thermo.cs2;
        (*dwds_array)[idx] = thermo.dwds;
        (*dwdB_array)[idx] = thermo.dwdB;
        (*dwdS_array)[idx] = thermo.dwdS;
        (*dwdQ_array)[idx] = thermo.dwdQ;

    }
    end = std::chrono::system_clock::now();
    cout << "EoS table generation: 100% done. Time elapsed: " << std::chrono::duration<double>(end-start).count() << "s" << endl;

    //==========================================================================
    //Write to HDF5 file
    //==========================================================================
    #ifdef SAVE_DOMAIN
    s_dataset.write(s_array->data(), H5::PredType::NATIVE_DOUBLE);
    rhoB_dataset.write(rhoB_array->data(), H5::PredType::NATIVE_DOUBLE);
    rhoS_dataset.write(rhoS_array->data(), H5::PredType::NATIVE_DOUBLE);
    rhoQ_dataset.write(rhoQ_array->data(), H5::PredType::NATIVE_DOUBLE);
    #endif
    T_dataset.write(T_array->data(), H5::PredType::NATIVE_DOUBLE);
    muB_dataset.write(muB_array->data(), H5::PredType::NATIVE_DOUBLE);
    muS_dataset.write(muS_array->data(), H5::PredType::NATIVE_DOUBLE);
    muQ_dataset.write(muQ_array->data(), H5::PredType::NATIVE_DOUBLE);
    e_dataset.write(e_array->data(), H5::PredType::NATIVE_DOUBLE);
    p_dataset.write(p_array->data(), H5::PredType::NATIVE_DOUBLE);
    cs2_dataset.write(cs2_array->data(), H5::PredType::NATIVE_DOUBLE);
    dwds_dataset.write(dwds_array->data(), H5::PredType::NATIVE_DOUBLE);
    dwdB_dataset.write(dwdB_array->data(), H5::PredType::NATIVE_DOUBLE);
    dwdS_dataset.write(dwdS_array->data(), H5::PredType::NATIVE_DOUBLE);
    dwdQ_dataset.write(dwdQ_array->data(), H5::PredType::NATIVE_DOUBLE);

    //==========================================================================
    //Close HDF5 file and release memory
    //==========================================================================
    #ifdef SAVE_DOMAIN
    delete s_array, rhoB_array, rhoS_array, rhoQ_array;
    #endif
    delete T_array, muB_array, muQ_array, muS_array,
           e_array, p_array, cs2_array,
           dwds_array, dwdB_array, dwdQ_array, dwdS_array; //Free memory
    }

    ~Inverter(){
      eos_h5_file.close();
      delete[] eos;
    }
};


int main(int argc, char** argv){

    if (check_args(argc, argv) != 0) return 1;
    //Set number of threads
    int nthreads = atoi(argv[1]);
    omp_set_num_threads(nthreads);

    Inverter inverter(argv[4], nthreads);
    inverter.init_ccake(argv[2], argv[3]);

    int Nsize = 5*4*4*4;
    auto start = std::chrono::system_clock::now();
    for (int is=-3; is<2; ++is){
      for (int iB=-3; iB<1; ++iB){
        for (int iS=-3; iS<1; ++iS){
          for (int iQ=-3; iQ<1; ++iQ){
            inverter.set_range({is, iB, iS, iQ});
            std::string tableName = "table_" + std::to_string(is) + "_"
                                             + std::to_string(iB) + "_"
                                             + std::to_string(iS) + "_"
                                             + std::to_string(iQ);
            int idx = (is+3)*4*4*4 + (iB+3)*4*4 + (iS+3)*4 + iQ+3;
            inverter.create_table(tableName);
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;
            std::cout << "EoS table generation: " << idx << " / " << Nsize << " = " << (double)idx/Nsize*100. << "% done. Elapsed time: " << elapsed_seconds.count() << "s";
            std::cout << " Estimated time remaining: " << elapsed_seconds.count()/(double)idx*(Nsize-idx) << "s" << std::endl;
          }
        }
      }
    }
    return 0;
}