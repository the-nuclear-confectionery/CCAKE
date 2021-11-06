#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "interpolator4D.h"

void Interpolator4D::initialize( string filename )
{
  // import data from files (move this to I/O class eventually)
  load_data( filename );

  // set up all relevant quantities needed for evaluation
  construct_interpolant();
}


void Interpolator4D::load_data( string filename )
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
        if ( dim != 4 )
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

void Interpolator4D::construct_interpolant()
{
  // initialize grid ranges, etc.
  for ( auto & gridDirection : grid )
  {
    grid_mins.push_back( *min_element(gridDirection.begin(), gridDirection.end()) );
    grid_maxs.push_back( *max_element(gridDirection.begin(), gridDirection.end()) );
  }

  const int dim = 4;
  grid_points.resize(dim);
  for ( int iDim = 0; iDim < dim; iDim++ )
  {
    const int grid_n = grid_sizes[iDim];
    const double delta = (grid_mins[iDim] - grid_mins[iDim]) / (grid_n-1.0);
    grid_spacings.push_back( delta );
    for ( int iGrid = 0; iGrid < grid_n; iGrid++ )
      grid_points[iDim].push_back( grid_mins[iDim] + iGrid*delta );
  }
}


void Interpolator4D::evaluate( vector<double> & coordinates, vector<double> & results )
{
  //const int dim = coordinates.size();
  const int dim = 4;

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
    cout << "CHECK: " << ic << "   " << fracs[ic] << "   " << inds[ic] << endl;

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
  for (int iT = 0; iT < 1; iT++)
  {
    for (int imuB = 0; imuB < 1; imuB++)
    {
      for (int imuS = 0; imuS < 1; imuS++)
      {
        for (int imuQ = 0; imuQ < 1; imuQ++)
        {
          double weight = fracs[0]*fracs[1]*fracs[2]*fracs[3]; // fractional weights
          auto & cell = fields[ indexer( inds[0] + iT,   inds[1] + imuB,
                                         inds[2] + imuS, inds[3] + imuQ ) ];
          // 
          for ( int iField = 0; iField < nFields; iField++ )
            results[iField] += weight * cell[iField];

          // invert each weight to get opposite cell's contribution
          fracs[3] = 1.0 - fracs[3];
        }
        fracs[2] = 1.0 - fracs[2];
      }
      fracs[1] = 1.0 - fracs[1];
    }
    fracs[0] = 1.0 - fracs[0];
  }
}