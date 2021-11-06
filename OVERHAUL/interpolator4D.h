#ifndef INTERPOLATOR4D_H
#define INTERPOLATOR4D_H

//#include <map>
#include <string>
#include <vector>

using namespace std;

class Interpolator4D
{
  public:
    Interpolator4D(){}
    ~Interpolator4D(){}

    //void set_columns( vector<string> & columns_in );
    void initialize( string filename );

    void evaluate( vector<double> & coordinates, vector<double> & results );

  private:
    vector<int> grid_sizes;
    vector<double> grid_spacings, grid_mins, grid_maxs;

    /*std::map<string,int> columns;
    first['a']=10;
    first['b']=30;*/
    vector<vector<double> > grid, grid_points, fields;

    inline size_t indexer( vector<int> & indices )
    {
      return ( ( indices[0] * grid_sizes[1] + indices[1] )
                * grid_sizes[2] + indices[2] )
                * grid_sizes[3] + indices[3];
    }

    inline size_t indexer( int & i0, int i1, int i2, int i3 )
    {
      return ( ( i0 * grid_sizes[1] + i1 ) * grid_sizes[2] + i2 ) * grid_sizes[3] + i3;
    }

    void load_data( string filename );
    void construct_interpolant();

};



#endif