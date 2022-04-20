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
  void initialize( double it0, int ntot, double h,
                   vector<Particle> * particlesPtr_in,
                   double dtsave, int & numpart );
  int triToSum( Vector<int,2> dael, Vector<int,2> size );
  //integers
  int _n;
  int n() { return _n; }
  int start         = 0;
  int end           = 0;
  int fnum          = 0;
  int qmf           = 0;
  //if==1 quantum mechanicanical corrections to the flow or added
  //if==0 no corrections are included
  //int number_part   = 0;
  //int rk2           = 0; 
  int steps         = 0;
  int gtyp          = 0;
  int first         = 1;  //HARDCODE
  int average       = 0;
  int lowT          = 0;
  int etaconst      = 0;
  int fcount        = 0;
  int cevent        = 0;
  int range         = 0; //range is number of boxes from left to right extra
  int Size          = 0;
  int cfon          = 1;  // HARDCODE
  int visc          = 0; 
  // visc=0 for ideal
  // visc=1 for bulk
  // visc=2 for shear
  // visc=3 for bulk+shear
  // visc=4 for BSQ+bulk+shear
  
  //doubles
  double _h         = 0;
  double gd2        = 0;
  double t0         = 0;
  double t          = 0;
  double dt         = 0;
  double factor     = 0;
  /*double dEz        = 0;
  double E          = 0;
  double Ez         = 0;
  double E0         = 0;
  double Etot       = 0;
  double S          = 0;
  double S0         = 0;
  double Eloss      = 0;
  double Esubb      = 0;
  double Btotal     = 0;
  double Stotal     = 0;
  double Qtotal     = 0;
  double Btotal0    = 0;
  double Stotal0    = 0;
  double Qtotal0    = 0;*/
  /*double bvf        = 0;
  double svf        = 0; 
  double zwidth     = 0;
  double sTc        = 0;
  double zTc        = 0;*/
  double E1         = 0;
  double E2         = 0;
  double step       = 0;
  double efcheck    = 0;
  double sfcheck    = 0;

  //vectors of int
  vector<int> list;
  vector<int> lead;
  vector<int> link;
  Vector<int,2> size;
  
  //vectors of vectors
  vector< Vector<int,2> > dael;

  //strings
  string eos_s       = "";
  string eost        = "";
  string eos_p       = "";
  string ebe_folder  = "";
  //vector of strings
  vector<string> filenames;

  //vector of pointers
  vector<Particle> * particlesPtr;

  static constexpr double tend   = 50.02;

  //private:

  static constexpr int VERBOSE  = 5;
  static constexpr double e0    = 1.0;
  static constexpr double q     = 1.0;
  static constexpr double dTemp = 0.00001;

  //vectors of doubles
  Vector<double,2> min;
  Vector<double,2> max;
  Vector<double,2> uni;

  double gradPressure_weight(vector<Particle> & particles, const int a, const int b)
  {
    const auto & pa = particles[a];
    const auto & pb = particles[b];

    double alpha_q = 1.0;
    double v_signal_q = sqrt(1.0/3.0);

    double innerp = inner( pa.r - pb.r, pa.qmom - pb.qmom );
    double innerr = inner( pa.r - pb.r, pa.r    - pb.r    );
    innerp = 2.0*alpha_q*v_signal_q
            / ( pa.sigma/pa.gamma + pb.sigma/pb.gamma )
            / sqrt(innerr) * innerp;

    if ( innerp > 0.0 || a==b ) innerp=0.0;


    return pb.sigmaweight*pa.sigma
          *(    pb.eosPtr->p() / (pb.sigma*pb.sigma)
              + pa.eosPtr->p() / (pa.sigma*pb.sigma)
              - innerp );
  }

};

#endif
