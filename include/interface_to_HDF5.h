#ifndef INTERFACE_TO_HDF5_H
#define INTERFACE_TO_HDF5_H

#include <iomanip>
//#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

using std::string;
using std::vector;

class interface_to_HDF5
{

  private:
    H5File file;
    string GROUPEVENT_NAME = "/Event";

  public:

    interface_to_HDF5(){}
    ~interface_to_HDF5(){}

    //--------------------------------------------------------------------------
    void output_double_attribute(Group & group, double value, string name)
    {
      hsize_t dims[1]; dims[0] = 1;
      DataSpace dspace(1, dims);
      Attribute att = group.createAttribute(name.c_str(), PredType::NATIVE_DOUBLE, dspace );
      att.write(PredType::NATIVE_DOUBLE, &value);

      return;
    }


    //--------------------------------------------------------------------------
    void initialize( const std::string & filename,
                     const vector<double> & parameters,
                     const vector<string> & parameter_names )
    {
      try
      {
        Exception::dontPrint();

        const H5std_string FILE_NAME(filename);
        GROUPEVENT_NAME = "/Event";

        file = H5File(FILE_NAME, H5F_ACC_TRUNC);
        Group groupEvent(file.createGroup(GROUPEVENT_NAME.c_str()));

        for (int i = 0; i < parameters.size(); i++)
          output_double_attribute( groupEvent, parameters[i],
                                               parameter_names[i] );

      }

      // catch any errors
      catch(FileIException error)
      {
        error.printErrorStack();
        abort();
      }

      catch(DataSetIException error)
      {
        error.printErrorStack();
        abort();
      }

      catch(DataSpaceIException error)
      {
        error.printErrorStack();
        abort();
      }

      return;
    }



    //------------------------------------------------------------------------------
    string get_zero_padded_int( int i, int width )
    {
      std::ostringstream ss;
        ss << std::setw( width ) << std::setfill( '0' ) << i;
        return ss.str();
    }



    void set_units(DataSet & ds, const std::string & units)
    {
      StrType str_type(PredType::C_S1, H5T_VARIABLE);
      DataSpace att_space(H5S_SCALAR);
      Attribute att = ds.createAttribute( "Units", str_type, att_space );
      att.write( str_type, units );
      return;
    }


    //------------------------------------------------------------------------------
    void output_dataset( const vector<string> & dataset_names,
                         const vector<string> & dataset_units,
                         const vector<vector<double> > & v, const int width,
                         const int timestepindex, const vector<int> & labels,
                         const vector<double> & parameters,
                         const vector<string> & parameter_names )
    {
      const int Nfields = v.size();
      const int length = v[0].size();

      // double data[Nfields][length];
      // for (int i = 0; i < Nfields; i++)
      // for (int j = 0; j < length;  j++)	// skip x and y coordinates
      //   data[i][j] = v[i][j];


      string FRAME_NAME = GROUPEVENT_NAME + "/Frame_"
                + get_zero_padded_int( timestepindex, width );

      Group groupFrame(file.createGroup(FRAME_NAME.c_str()));

      //output_double_attribute( groupFrame, time, "Time" );
      for (int i = 0; i < parameters.size(); i++)
        output_double_attribute( groupFrame, parameters[i],
                                             parameter_names[i] );

      // export all datasets
      for (int iDS = 0; iDS < dataset_names.size(); iDS++)
      {
        hsize_t dims[1];
        dims[0] = length;
        DataSpace dataspace(1, dims);

        string DATASET_NAME = FRAME_NAME + "/" + dataset_names[iDS];
        DataSet dataset = groupFrame.createDataSet( DATASET_NAME.c_str(),
                              PredType::NATIVE_DOUBLE, dataspace);

        // specify units
        set_units( dataset, dataset_units[iDS] );

        dataset.write(v[iDS].data(), PredType::NATIVE_DOUBLE);
      }

      // export any particle labels as indices (e.g., current eos names)
      {
        hsize_t dims[1];
        dims[0] = length;
        DataSpace dataspace(1, dims);

        string DATASET_NAME = FRAME_NAME + "/labels";
        DataSet dataset = groupFrame.createDataSet( DATASET_NAME.c_str(),
                              PredType::NATIVE_INT, dataspace);

        dataset.write(labels.data(), PredType::NATIVE_INT);
      }

//if (1) exit(-1);

      return;
    }


//    // SAVE THIS FOR NOW
//    // a typedef for our managed H5File pointer
//    typedef std::shared_ptr<H5::H5File> H5FilePtr;
//
//    // create or open a file
//    H5FilePtr create_or_open(const std::string& fname)
//    {
//        H5::Exception::dontPrint();
//        H5::H5File* file = 0;
//
//        try {
//            file = new H5::H5File(fname.c_str(), H5F_ACC_RDWR);
//        } catch(const H5::FileIException&) {
//            file = new H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
//        }
//
//        return H5FilePtr(file);
//    }

};

#endif
