#include <algorithm>
#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../include/interpolatorND.h"

////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::initialize( string filename, double fill_value_in )
{
  fill_value = fill_value_in;

  // need this to iterate systematically over hypercube vertices
  initialize_hypercube();

  // import data from files (move this to I/O class eventually)
  load_data( filename );

  // set up all relevant quantities needed for evaluation
  construct_interpolant();
}


////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::load_data( string filename )
{
  if ( using_HDF )
    load_data_from_HDF( filename );
  else
    load_data_from_dat( filename );
}

////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::load_data_from_HDF( string filename )
{
  H5File f( filename, H5F_ACC_RDONLY );

  // Reading in grid dimensions first
  hsize_t gridsizes_dims[1];

  DataSet ds = f.openDataSet( "/dimensions" );
  DataSpace dspace = ds.getSpace().getSimpleExtentDims(gridsizes_dims, NULL);

  const int dimension_of_grid = gridsizes_dims[0];
  cout << "Dims: " << dimension_of_grid << endl;

  int grid_dimensions[ dimension_of_grid ];
  ds.read(grid_dimensions, H5::PredType::NATIVE_INT);

  cout << "Grid dimensions:";
  for (int iDim = 0; iDim < dimension_of_grid; iDim++)
    cout << " " << grid_dimensions[iDim];
  cout << endl;

  // must set grid_sizes for below
  grid_sizes.assign(grid_dimensions, grid_dimensions + dimension_of_grid);

  // Reading in data grid next
  hsize_t data_dims[2];

  DataSet ds2 = f.openDataSet( "/data" );
  DataSpace dspace2 = ds2.getSpace().getSimpleExtentDims(data_dims, NULL);

  const size_t nRows = data_dims[0];
  const size_t nCols = data_dims[1];
  cout << "Dims: " << nRows << " " << nCols << endl;

  cout << "Creating data array" << endl;

  double * data = new double [ nRows * nCols ];

  cout << "Created data array" << endl;

  // read in the data
  ds2.read(data, H5::PredType::NATIVE_DOUBLE);

  // define convenient lambda to index data as 1D array
  auto H5_indexer = [nCols](size_t i, size_t j) { return i*nCols+j; };

  //============================================================================
  // FINALLY, need to set interpolator quantities as in load_data_from_dat()

  // check dimensions
  if ( dimension_of_grid != D )
  {
    std::cerr << "Your equation of state table has the wrong dimension!" << std::endl;
    abort();
  }

  // "grid" holds coordinates of grid points
  grid.resize(dimension_of_grid);
  for ( int iCol = 0; iCol < dimension_of_grid; iCol++ )
  {
    grid[iCol].resize( nRows );
    for ( size_t iRow = 0; iRow < nRows; iRow++ )
      grid[iCol].at( iRow ) = data[ H5_indexer( iRow, iCol ) ];
  }

  // "fields" holds all quantities to be interpolated over
  fields.resize( nRows );
  for ( size_t iRow = 0; iRow < nRows; iRow++ )
  {
    // quantities are whatever columns are left after getting the grid points
    fields[iRow].resize( nCols - dimension_of_grid );
    for ( int iCol = dimension_of_grid; iCol < nCols; iCol++ )
      fields[iRow].at( iCol - dimension_of_grid )
        = data[ H5_indexer( iRow, iCol ) ];
  }

  delete [] data;
}




////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::load_data_from_dat( string filename )
{
  ifstream infile( filename.c_str() );
  if (infile.is_open())
  {
    string line;

    int n_header_line = 1;  // assume one header line containing grid sizes
    int count = 0;
    int dim = 0;            // use this to enable more general ND interpolator
    while ( getline (infile, line) )
    {
      istringstream iss(line);
      if ( count++ < n_header_line )
      {
        int tmp = 0;
        while (iss >> tmp) grid_sizes.push_back(tmp);
        dim = grid_sizes.size();
        if ( dim != D )
        {
          std::cerr << "Your equation of state table has the wrong dimension!" << std::endl;
          abort();
        }
        grid.resize(dim);
      }
      else
      {
        ///////////////////////////////
        // set grid coordinates first
        double tmp = 0.0;
        for ( int iDim = 0; iDim < dim; iDim++ )
        {
          iss >> tmp;
          grid[iDim].push_back( tmp );
        }

        ///////////////////////////////
        // the rest of the columns are fields to interpolate
        vector<double> field;
        while (iss >> tmp) field.push_back(tmp);
        fields.push_back( field );
      }
    }

//    cout << "Test grid integrity" << endl;
//    bool crash = false;
//    size_t tmpsize = fields[0].size();
//    for (auto&field:fields)
//      if (field.size()!=tmpsize)
//      {
//        crash = true;
//        cout << field.size() << "=!=" << tmpsize << endl;
//      }
//
//    if (crash) cout << "Failed!" << endl; else cout << "Succeeded!" << endl;
//    cout << "tmpsize = " << tmpsize << endl;
//    abort();

    infile.close();
  }
  else
  {
    std::cerr << "File " << filename << " could not be opened!\n";
    abort();
  }
}


////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::construct_interpolant()
{
  grid_mins.clear();
  grid_maxs.clear();
  grid_points.clear();
  grid_spacings.clear();

  // initialize grid ranges, etc.
  for ( auto & gridDirection : grid )
  {
    grid_mins.push_back( *min_element(gridDirection.begin(), gridDirection.end()) );
    grid_maxs.push_back( *max_element(gridDirection.begin(), gridDirection.end()) );
  }

  // using grid ranges and sizes, re-construct grid lattices themselves
  const int dim = D;
  grid_points.resize(dim);
  for ( int iDim = 0; iDim < dim; iDim++ )
  {
    const int grid_n = grid_sizes[iDim];
    const double delta = (grid_maxs[iDim] - grid_mins[iDim]) / (grid_n-1.0);
    grid_spacings.push_back( delta );
    for ( int iGrid = 0; iGrid < grid_n; iGrid++ )
      grid_points[iDim].push_back( grid_mins[iDim] + iGrid*delta );
  }

}

////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::initialize_hypercube()
{
  const int dim = D;
  n_points_in_hypercube = 1;
  for ( int iDim = 0; iDim < dim; iDim++ ) n_points_in_hypercube *= 2;

  for ( int hypercube_index = 0; hypercube_index < n_points_in_hypercube; hypercube_index++ )
  {
    bitset<D> bs( hypercube_index );

    vector<int> digits;
    for (int iDigit = bs.size()-1; iDigit >= 0; iDigit--)
      digits.push_back( bs[iDigit] );

    hypercube_indices.push_back( digits );
  }
}

////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::evaluate( const vector<double> & coordinates,
                                  vector<double> & results )
{
  // if coordinates are nan's, don't bother evaluating!
  for (auto & c: coordinates)
    if ( isnan(c) )
    {
      results = vector<double>(fields.front().size(), fill_value);
      return;
    }

  const int dim = D;

  bool out_of_range = false;

  //////////////////////////////////////////
  // locate coordinate in grid
  vector<int> inds(dim);     // integral coordinate location for containing hypercube
  vector<double> fracs(dim); // fractional coordinate location in containing hypercube
  for ( int ic = 0; ic < dim; ic++ )
  {
    double index = 0.0;      // holds integer part of the coordinate location in grid
    fracs[ic] = 1.0 - modf( (coordinates[ic] - grid_mins[ic])
                            / grid_spacings[ic], &index );
    inds[ic] = static_cast<int>( index );
//    cout << "CHECK: " << ic << "   " << fracs[ic] << "   " << inds[ic] << "   " << index
//        << "   " << coordinates[ic] << "   " << (coordinates[ic] - grid_mins[ic])
//                                                / grid_spacings[ic]
//        << "   " << 1.0 - fracs[ic] + index << endl;

    // handle special situation where queried point is at the grid maximum
    if ( inds[ic] + 1 == grid_sizes[ic] )
    {
      // reverse interpolation weights along this dimension
      fracs[ic] = 1.0 - fracs[ic];
      inds[ic]--;
    }

    // check if this index is out of range
    if ( inds[ic] < 0 || inds[ic] >= grid_sizes[ic] )
      out_of_range = true;
  }

  //////////////////////////////////////////
  // compute linear interpolant
  // all fields interpolated at once
  const int nFields = fields.front().size();
  if ( out_of_range )
  {
    results = vector<double>(nFields, fill_value);
    return;
  }
  else
    results = vector<double>(nFields, 0.0);

//  cout << "coords:";
//  for (auto & c: coordinates)
//    cout << " " << c;
//  cout << endl;

  // loop over hypercube indices
  for ( auto & hypercube_index : hypercube_indices )
  {
    double weight = 1.0;
    for (int iDim = 0; iDim < dim; iDim++)
      weight *= (fracs[iDim] + hypercube_index[iDim]
                  - 2.0*hypercube_index[iDim]*fracs[iDim]);

    vector<int> hypercube_inds = inds;
    for (int iDim = 0; iDim < dim; iDim++)
      hypercube_inds[iDim] += hypercube_index[iDim];

//    cout << "NEED TO FIX ISSUES BEFORE USING THIS" << endl;
//    if (1) abort();

    //cout << "inds:";
    //for (auto&is:inds) cout << " " << is;
    //cout << endl;
//    cout << "hypercube_inds:";
//    for (auto&his:hypercube_inds) cout << " " << his;
//    cout << endl;
//    cout << "fields.size() = " << fields.size() << endl;
//    cout << "indexer( hypercube_inds ) = " << indexer( hypercube_inds ) << endl;
    auto & cell = fields[ indexer( hypercube_inds ) ];

    //cout << "results.size() = " << results.size() << endl;
    //cout << "cell.size() = " << cell.size() << endl;
//    if (cell.size() != 17)
//    {
//      //cout << "Test grid integrity" << endl;
//      bool crash = false;
//      size_t tmpsize = fields[0].size();
//      for (auto&field:fields)
//        if (field.size()!=tmpsize)
//        {
//          crash = true;
//          cout << field.size() << "=!=" << tmpsize << endl;
//        }
//
//      //if (crash) cout << "Failed!" << endl; else cout << "Succeeded!" << endl;
//      //cout << "tmpsize = " << tmpsize << endl;
//      //abort();
//    }
    // 
    for ( int iField = 0; iField < nFields; iField++ )
    {
//      if (iField==6)
//      {
//        for (auto&his:hypercube_inds) cout << " " << his;
//        cout << "   " << cell[iField] << endl;
//      }
//      cout << "Check nodes: " << iField << "   " << weight << "   " << cell[iField] << endl;
      results[iField] += weight * cell[iField];
    }
  }
  return;
}




////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::evaluate(
      const vector<double> & coordinates, vector<double> & results,
      const vector<string> & fields_to_interpolate )
{
  // if coordinates are nan's, don't bother evaluating!
  for (auto & c: coordinates)
    if ( isnan(c) )
    {
      results = vector<double>(fields_to_interpolate.size(), fill_value);
      return;
    }

  const int dim = D;

//  cout << "-----------------------------------------" << endl;
//  for (auto&coord:coordinates) cout << " " << coord;
//  cout << endl;
//  cout << "-----------------------------------------" << endl;


  bool out_of_range = false;

  //////////////////////////////////////////
  // locate coordinate in grid
  vector<int> inds(dim);     // integral coordinate location for containing hypercube
  vector<double> fracs(dim); // fractional coordinate location in containing hypercube
  for ( int ic = 0; ic < dim; ic++ )
  {
    double index = 0.0;      // holds integer part of the coordinate location in grid
    fracs[ic] = 1.0 - modf( (coordinates[ic] - grid_mins[ic])
                            / grid_spacings[ic], &index );
    inds[ic] = static_cast<int>( index );
//    cout << "CHECK: " << ic << "   " << fracs[ic] << "   " << inds[ic] << "   " << index
//        << "   " << coordinates[ic] << "   " << (coordinates[ic] - grid_mins[ic])
//                                                / grid_spacings[ic] << endl;

    // handle special situation where queried point is at the grid maximum
    if ( inds[ic] + 1 == grid_sizes[ic] )
    {
      // reverse interpolation weights along this dimension
      fracs[ic] = 1.0 - fracs[ic];
      inds[ic]--;
    }

    // check if this index is out of range
    if ( inds[ic] < 0 || inds[ic] >= grid_sizes[ic] )
      out_of_range = true;
  }

  //////////////////////////////////////////
  // compute linear interpolant
  // all fields interpolated at once
  const int nFields = fields_to_interpolate.size();
  if ( out_of_range )
  {
    results = vector<double>(nFields, fill_value);
    return;
  }
  else
    results = vector<double>(nFields, 0.0);

//  cout << "coords:";
//  for (auto & c: coordinates)
//    cout << " " << c;
//  cout << endl;

  // set field indices here (map is slow)
  vector<int> field_indices(nFields, -1);
  for ( int iField = 0; iField < nFields; iField++ )
    field_indices[iField] = field_names[fields_to_interpolate[iField]];

  // loop over hypercube indices
  for ( auto & hypercube_index : hypercube_indices )
  {
    double weight = 1.0;
    for (int iDim = 0; iDim < dim; iDim++)
      weight *= (fracs[iDim] + hypercube_index[iDim]
                  - 2.0*hypercube_index[iDim]*fracs[iDim]);

    vector<int> hypercube_inds = inds;
    for (int iDim = 0; iDim < dim; iDim++)
      hypercube_inds[iDim] += hypercube_index[iDim];

//    cout << "inds:";
//    for (auto&is:inds) cout << " " << is;
//    cout << endl;
//    cout << "hypercube_inds:";
//    for (auto&his:hypercube_inds) cout << " " << his;
//    cout << endl;
//    cout << "fields.size() = " << fields.size() << endl;
//    cout << "indexer( hypercube_inds ) = " << indexer( hypercube_inds ) << endl;

    auto & cell = fields[ indexer( hypercube_inds ) ];

    // 
    for ( int iField = 0; iField < nFields; iField++ )
    {
//      if (iField==6)
//      {
//        for (auto&his:hypercube_inds) cout << " " << his;
//        cout << "   " << cell[field_names[fields_to_interpolate[iField]]] << endl;
//      }
      results[iField] += weight * cell[ field_indices[iField] ];
    }
  }
  return;
}




////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::rescale( const string & column_to_rescale,
                                 const string & column_to_rescale_by,
                                 int power_of_rescaling )
{
  cout << "Rescaling " << column_to_rescale << " by "
      << column_to_rescale_by << " to the power of " << power_of_rescaling << endl;
  int column_index_to_rescale    = field_names[column_to_rescale];
  int column_index_to_rescale_by = grid_names[column_to_rescale_by];
  auto & grid_column             = grid[column_index_to_rescale_by];

  const size_t nCells = fields.size();
  for (size_t iCell = 0; iCell < nCells; iCell++)
    fields[iCell][column_index_to_rescale]
      *= pow( grid_column[iCell], power_of_rescaling );

  return;
}

////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::rescale_axis( const string & column_to_rescale,
                                      double overall_factor )
{
  cout << "Rescaling " << column_to_rescale << " by factor of " << overall_factor << endl;
  int column_index_to_rescale = grid_names[column_to_rescale];

  auto & axis = grid[column_index_to_rescale];
  for ( auto & pt : axis ) pt *= overall_factor;

  // reset all quantities appropriately (not efficient for each axis individually)
  construct_interpolant();

  return;
}

////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::rescale_axes( double overall_factor )
{
  // rescale grid by factor
  for ( auto & axis : grid )
  for ( auto & pt : axis )
    pt *= overall_factor;

  // reset all quantities appropriately
  construct_interpolant();

  return;
}

////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::rescale_field( const string & column_to_rescale,
                                       double overall_factor )
{
  int column_index_to_rescale = field_names[column_to_rescale];
  for ( auto & cell : fields ) cell[column_index_to_rescale] *= overall_factor;
  return;
}


////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::apply_function_to_field(
        const string & column_to_modify, std::function<double(double)> f )
{
  int column_index_to_modify = field_names[column_to_modify];
  for ( auto & cell : fields )
    cell[column_index_to_modify] = f( cell[column_index_to_modify] );
  return;
}
