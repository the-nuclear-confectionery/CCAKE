#ifndef LINKLIST_H
#define LINKLIST_H

class LinkList
{
public:
  LinkList();
  ~LinkList();

  static constexpr double tend=50.02;

  double _h;
  int _n;
  int start, end, fnum;
  vector<double> sFO, Tfluc; //entropy at freezeout
  int qmf; //if==1 quantum mechanicanical corrections to the flow or added
           //if==0 no corrections are included

  double kernel( Vector<double,2> kin );

  double knorm, knorm2, kgrad, kgrad2;
  int number_part;
  double t0;
  double t, dt;
  double factor;
  int frzc;
  double tau, taup, taupp;
  int rk2;
  double gd2;

  int gtyp;

  int cfon;
  vector<int> list;
  //int *lead;
  //int *link;
  vector<int> lead;
  vector<int> link;
  int cf;
  int visc; // visc=0 for ideal
            // visc=1 for bulk
            // visc=2 for shear
            // visc=3 for bulk+shear
  double efcheck, sfcheck;
  int steps;

  int first;
  int average;
  int lowT;

  double *divTtemp, *gsub, *bulksub, *swsub, *shear33sub, *tlist;
  double avgetasig; // possibly not needed?
  Matrix<double,3,3> *shearsub;
  Vector<double,2> *divT,*rsub;
  Vector<double,2> *uout;
  double wfz,cs2;

  Vector<int,2> size;

  //Vector<int,2> *dael;
  vector< Vector<int,2> > dael;

  double dEz, E, Ez, E0, Etot, S, S0, Eloss, Esubb;
  double Btotal, Stotal, Qtotal;
  double Btotal0, Stotal0, Qtotal0;

  double E1, E2;


  int etaconst;
  double bvf, svf, zwidth, sTc, zTc;

  void setup( double it0, int ntot, double h,
              vector<Particle> & particles, double dtsave, int & numpart );

  string eos_s,eos_p;
  vector<string> filenames;

  string ebe_folder;

  int fcount,cevent;
  string eost;

  void initialize();

  int triToSum( Vector<int,2> dael, Vector<int,2> size );

  void destroy();

  int n() { return _n; }

  void etas_set();
  void bsqsv_set();
  void endEV();
  void setshear();
  void prints(); //possibly not needed?


private:

  static constexpr int VERBOSE  = 5;
  static constexpr double e0    = 1.0;
  static constexpr double q     = 1.0;
  static constexpr double dTemp = 0.00001;

  int range;                //range is number of boxes from left to right extra
  int Size;
  double step;

  Vector<double,2> min;
  Vector<double,2> max;
  Vector<double,2> uni;

  double gradPressure_weight(int a, int b)
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
          *(    pb.eos.p() / (pb.sigma*pb.sigma)
              + pa.eos.p() / (pa.sigma*pb.sigma)
              - innerp );
  }

  Vector<double,2> gradKernel( Vector<double,2> a, bool verbose = false );

}

#endif