#include <algorithm>
#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "interpolatorND.h"

template <int D>
void InterpolatorND<D>::initialize( string filename )
{
  // import data from files (move this to I/O class eventually)
  load_data( filename );

  // set up all relevant quantities needed for evaluation
  construct_interpolant();
}


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
        //size_t grid_length = 1;
        //for ( auto & size : grid_sizes ) grid_length *= size;
        grid.resize(dim);
      }
      else
      {
        ///////////////////////////////
        // set grid coordinates first
        //vector<double> gridPoint(dim);
        //for ( auto & coord : gridPoint ) iss >> coord;
        //grid.push_back( gridPoint );
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

    infile.close();
  }
}

template <int D>
void InterpolatorND<D>::construct_interpolant()
{
  // initialize grid ranges, etc.
  for ( auto & gridDirection : grid )
  {
    grid_mins.push_back( *min_element(gridDirection.begin(), gridDirection.end()) );
    grid_maxs.push_back( *max_element(gridDirection.begin(), gridDirection.end()) );
  }

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


template <int D>
void InterpolatorND<D>::evaluate( vector<double> & coordinates, vector<double> & results )
{
  //const int dim = coordinates.size();
  const int dim = D;

  //////////////////////////////////////////
  // locate coordinate in grid
  vector<int> inds(dim);     // integral coordinate location for containing hypercube
  vector<double> fracs(dim); // fractional coordinate location in containing hypercube
  for ( int ic = 0; ic < dim; ic++ )
  {
    double index = 0.0; // holds integer part of the coordinate location in grid
    fracs[ic] = 1.0 - modf( (coordinates[ic] - grid_mins[ic])
                            / grid_spacings[ic], &index );
    inds[ic] = static_cast<int>( index );
    cout << "CHECK: " << ic << "   " << fracs[ic] << "   " << inds[ic] << "   " << index
        << "   " << (coordinates[ic] - grid_mins[ic]) / grid_spacings[ic] << endl;

  }

  //////////////////////////////////////////
  // compute linear interpolant (assuming 4D thermodynamics for simplicity)
  // all fields interpolated at once
  if (false)
  {
    std::cout << "You still need to check the order of the loops!" << endl;
    std::cerr << "You still need to check the order of the loops!" << endl;
    std::cout << "Also, you should invert fields if you haven't done so already" << std::endl;
    exit(1);
  }
  const int nFields = fields.front().size();
  results = vector<double>(nFields, 0.0);
  cout << "nFields = " << nFields << endl;

  // loop over hypercube indices
  for ( auto & hypercube_index : hypercube_indices )
  {
    double weight = 1.0;
    for (int iDim = 0; iDim < dim; iDim++)
      weight *= (fracs[iDim] + hypercube_index[iDim] - 2.0*hypercube_index[iDim]*fracs[iDim]);

    vector<int> hypercube_inds = inds;
    for (int iDim = 0; iDim < dim; iDim++)
      hypercube_inds[iDim] += hypercube_index[iDim];

    auto & cell = fields[ indexer( hypercube_inds ) ];

    // 
    for ( int iField = 0; iField < nFields; iField++ )
      results[iField] += weight * cell[iField];
  }
}