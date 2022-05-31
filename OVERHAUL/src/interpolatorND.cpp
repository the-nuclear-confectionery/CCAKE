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
void InterpolatorND<D>::initialize( string filename )
{
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
          exit(1);
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
  const int dim = D;

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
    cout << "CHECK: " << ic << "   " << fracs[ic] << "   " << inds[ic] << "   " << index
        << "   " << coordinates[ic] << "   " << (coordinates[ic] - grid_mins[ic])
                                                / grid_spacings[ic] << endl;
  }

  //////////////////////////////////////////
  // compute linear interpolant (assuming 4D thermodynamics for simplicity)
  // all fields interpolated at once
  const int nFields = fields.front().size();
  results = vector<double>(nFields, 0.0);
  cout << "nFields = " << nFields << endl;

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

    cout << "NEED TO FIX ISSUES BEFORE USING THIS" << endl;
    if (1) abort();

    cout << "inds:";
    for (auto&is:inds) cout << " " << is;
    cout << endl;
    cout << "hypercube_inds:";
    for (auto&his:hypercube_inds) cout << " " << his;
    cout << endl;
    cout << "fields.size() = " << fields.size() << endl;
    cout << "indexer( hypercube_inds ) = " << indexer( hypercube_inds ) << endl;
    auto & cell = fields[ indexer( hypercube_inds ) ];

    cout << "results.size() = " << results.size() << endl;
    cout << "cell.size() = " << cell.size() << endl;
    if (cell.size() != 17)
    {
      cout << "Test grid integrity" << endl;
      bool crash = false;
      size_t tmpsize = fields[0].size();
      for (auto&field:fields)
        if (field.size()!=tmpsize)
        {
          crash = true;
          cout << field.size() << "=!=" << tmpsize << endl;
        }

      if (crash) cout << "Failed!" << endl; else cout << "Succeeded!" << endl;
      cout << "tmpsize = " << tmpsize << endl;
      abort();
    }
    // 
    for ( int iField = 0; iField < nFields; iField++ )
      results[iField] += weight * cell[iField];
  }
}




////////////////////////////////////////////////////////////////////////////////
template <int D>
void InterpolatorND<D>::evaluate(
      const vector<double> & coordinates, vector<double> & results,
      const vector<string> & fields_to_interpolate )
{
  const int dim = D;

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
  }

  //////////////////////////////////////////
  // compute linear interpolant (assuming 4D thermodynamics for simplicity)
  // all fields interpolated at once
  const int nFields = fields_to_interpolate.size();
  results = vector<double>(nFields, 0.0);

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

    auto & cell = fields[ indexer( hypercube_inds ) ];

    // 
    for ( int iField = 0; iField < nFields; iField++ )
      results[iField] += weight * cell[field_names[fields_to_interpolate[iField]]];
  }
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