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
    int first         = 1;  //HARDCODE
    int range         = 0; //range is number of boxes from left to right extra
    int Size          = 0;
    
    //doubles
    double h         = 0;

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

    int n() { return n_particles; }

  private:

    static constexpr int VERBOSE  = 5;

    //vectors of doubles
    Vector<double,2> uni, min, max;

};

#endif
