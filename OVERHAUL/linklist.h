#ifndef LINKLIST_H
#define LINKLIST_H

#include "eos.h"
#include "particle.h"

class LinkList
{
  public:
    LinkList();
    ~LinkList(){}

    void reset();
    void initialize( vector<Particle> * particlesPtr_in, double h_in );
    inline int indexer( Vector<int,2> indices, Vector<int,2> sizes )
                { return indices(0) + indices(1)*sizes(0); }

    //integers
    int n_particles;
    int n() { return n_particles; }
//    int start         = 0;
//    int end           = 0;
//    int fnum          = 0;
//    int steps         = 0;
//    int gtyp          = 0;
    int first         = 1;  //HARDCODE
//    int average       = 0;
//    int lowT          = 0;
//    int etaconst      = 0;
    int range         = 0; //range is number of boxes from left to right extra
    int Size          = 0;
//    int cfon          = 1;  // HARDCODE
//    int visc          = 0; 
    // visc=0 for ideal
    // visc=1 for bulk
    // visc=2 for shear
    // visc=3 for bulk+shear
    // visc=4 for BSQ+bulk+shear
    
    //doubles
    double h         = 0;
//    double t0         = 0;
//    double t          = 0;
//    double dt         = 0;
//    double factor     = 0;

//    double E1         = 0;
//    double E2         = 0;
//    double step       = 0;
//    double efcheck    = 0;
//    double sfcheck    = 0;

    //vectors of int
    vector<int> list;
    vector<int> lead;
    vector<int> link;
    Vector<int,2> size;
    
    //vectors of vectors
    vector< Vector<int,2> > dael;
    vector<vector<int> > all_neighbors;

    //vector of pointers
    vector<Particle> * particlesPtr;


  private:

    static constexpr int VERBOSE  = 5;
//    static constexpr double e0    = 1.0;
//    static constexpr double q     = 1.0;

    //vectors of doubles
    Vector<double,2> min;
    Vector<double,2> max;
    Vector<double,2> uni;

};

#endif
