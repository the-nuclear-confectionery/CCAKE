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
    void initialize_hypercube();

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
    void rescale( const string & column_to_rescale, const string & column_to_rescale_by,
                  int power_of_rescaling );
    void rescale_axis( const string & column_to_rescale, double overall_factor );
    void rescale_axes( double overall_factor );
    void rescale_field( const string & column_to_rescale, double overall_factor );
    

    // interpolates all fields at once
    void evaluate( const vector<double> & coordinates, vector<double> & results );

    // interpolates only specific fields
    //void evaluate( const vector<double> & coordinates, vector<double> & results,
    //               vector<int> & fields_to_interpolate );

    // interpolates only specific fields by name
    void evaluate( const vector<double> & coordinates, vector<double> & results,
                   const vector<string> & fields_to_interpolate );

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
//      cout << "-------" << endl;
//      cout << grid_sizes.size() << " " << indices.size() << endl;
//      cout << "grid_sizes:" << endl;
//      for (auto&gs:grid_sizes) cout << " " << gs;
//      cout << endl << "-------" << endl;
//      for (auto&is:indices) cout << " " << is;
//      cout << endl << "-------" << endl;
      size_t result = indices[0];
      for ( size_t ind = 1; ind < indices.size(); ind++ )
      {
//        cout << ind << " " << grid_sizes.size() << " " << indices.size() << endl;
//        cout << grid_sizes[ind] << endl;
//        cout << indices[ind] << endl;
//        cout << result << endl;
//        cout << "-------" << endl;
        result = result * grid_sizes[ind] + indices[ind];
      }

      return result;
    }

    void load_data( string filename );
    void construct_interpolant();

};

#include "../src/interpolatorND.cpp"

#endif