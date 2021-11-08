#ifndef InterpolatorND_H
#define InterpolatorND_H

#include <map>
#include <string>
#include <vector>

using namespace std;

template <int D>
class InterpolatorND
{
  public:
    InterpolatorND<D>(){}
    ~InterpolatorND<D>(){}

    void initialize( string filename );

    void set_grid_names( const vector<string> & names )
    {
      int iName = 0;
      for ( auto & name : names ) grid_names[name] = iName++;
    }

    void set_field_names( const vector<string> & names )
    {
      int iName = 0;
      for ( auto & name : names ) field_names[name] = iName++;
    }

    // need this in case units are not physical (e.g., file contains dimensionless ratios)
    // Example: rescale( "e", "T", 4, hbarc_MeVfm )
    void rescale( string & column_to_rescale, string & column_to_rescale_by, int power_of_rescaling );
    void rescale_axis( string & column_to_rescale, double overall_factor );
    void rescale_field( string & column_to_rescale, double overall_factor );
    

    // interpolates all fields at once
    void evaluate( const vector<double> & coordinates, vector<double> & results );

    // interpolates only specific fields
    //void evaluate( const vector<double> & coordinates, vector<double> & results,
    //               vector<int> & fields_to_interpolate );

    // interpolates only specific fields by name
    void evaluate( const vector<double> & coordinates, vector<double> & results,
                   vector<string> & fields_to_interpolate );

    vector<double> get_grid_minima() { return grid_mins; }
    vector<double> get_grid_maxima() { return grid_maxs; }

  private:
    map<string, int> grid_names;
    map<string, int> field_names;

    int n_points_in_hypercube;
    vector<int> grid_sizes;
    vector<double> grid_spacings, grid_mins, grid_maxs;

    vector<vector<int> > hypercube_indices;
    vector<vector<double> > grid, grid_points, fields;

    inline size_t indexer( vector<int> & indices )
    {
      size_t result = indices[0];
      for ( size_t ind = 1; ind < indices.size(); ind++ )
        result = result * grid_sizes[ind] + indices[ind];
      return result;
    }

    void load_data( string filename );
    void construct_interpolant();

};

#include "interpolatorND.cpp"

#endif