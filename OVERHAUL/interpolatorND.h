#ifndef InterpolatorND_H
#define InterpolatorND_H

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

    void evaluate( vector<double> & coordinates, vector<double> & results );

  private:
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



#endif