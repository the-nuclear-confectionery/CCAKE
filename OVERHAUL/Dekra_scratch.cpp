
using namespace constants.h;

//=============================================================================
void read_in_file(string filename, vector<double> & event)
{
event.clear();
ifstream infile( filename.c_str() );
if (infile.is_open())
{
string line;
double phi, eta; 
while ( getline (infile, line) )
{
istringstream iss(line);
iss >> phi >> eta;
// impose cuts on eta
if ( abs(eta) > 2.4 ) continue;
event.push_back(phi + random_phi);
event.push_back(eta);
}
}
infile.close();
return;
}
//=============================================================================
void input_output::read_in_initial_conditions()
{
void readICs_iccing(string filename, int &_Ntable3, double factor, double const & efcheck,int &
 numpart)
{
//const double hbarC = 0.1973;
cout << "Reading in ICCING initial conditions!" << endl;
//string filename;
filename = ifolder+"ic0.dat";
ifstream infile( filename.c_str() );
cout << "Initial conditions file: " << filename << endl;
if (infile.is_open())
{
string line;
vector<double> xsub,ysub,esub,rBsub,rSsub,rQsub;
int count=0;
 while ( getline (infile, line) )
{
  istringstream iss(line);
  if count<n_header_lines
    { 
     iss >> xsub >> ysub >> esub >> rBsub >> rSsub >> rQsub;
    }
    else 
    {
    iss >> xsub >> ysub >> esub >> rBsub >> rSsub >> rQsub;  
    }

count ++;

if (!infile.is_open())
    {
        cout << "Can't open " << filename << endl;
        exit(1);
    }
  }
}

//double stepx,stepy;
//stringstream s;
//s << gx[1];
//s >> stepx;
stringstream s1;
    s1 << gx[2];
    s1 >> stepy;

cout << "dx=dy=" << stepx << " " << stepy << endl;
while (getline(input,line)) {
        std::vector<double> y (6,0) ;

        std::vector<std::string> x = split(line, ' ');


        for(int j=0; j<6; j++)
        {
            stringstream ss;
            ss << x[j];
            ss >> y[j];
//cout << "CHECK(" << __LINE__ << "): " << j << "   " << x[j] << "   " << y[j] << endl;
        }
}
}
}

//getline(input,line);
//std::vector<std::string> gx = split(line, ' ');
//==============================================================
// Add new check here to enforce freeze-out criterion before
// setting particle list size!!!

// do not scale by factor!!!
//if ((factor*y[2])>0.01)
if (y[2]>0.00301)
//if (y[2]>max(0.01,0.5*efcheck*hbarC)) //N.B. leave wiggle room around FO temp
{
// check if this file has finite energy density and zero charge densities
// (seems to cause problems in EOS)

const double eps_local = 1e-6;

/*if (   y[2] >= eps_local && y[3] < eps_local
	&& y[4] < eps_local  && y[5] < eps_local)
	{
		 y[3] = eps_local;
		 y[4] = eps_local;
		 y[5] = eps_local;
	}*/
xsub.push_back(y[0]);
ysub.push_back(y[1]);
esub.push_back(y[2]);
rBsub.push_back(y[3]);
rSsub.push_back(y[4]);
rQsub.push_back(y[5]);
}

input.close();

_Ntable3=xsub.size();
_p= new Particle<2>[_Ntable3];

cout << "After e-cutoff and freeze-out: size = " << _Ntable3 << endl;

int kk=_Ntable3;
numpart=0;	//number of frozen out particles

for(int j=0; j<_Ntable3; j++)
	{
const auto & pj = pj;
pj.r.x[0]=xsub[j];
pj.r.x[1]=ysub[j];
// pj.e_sub=EOS.e_out(factor*esub[j]);
pj.e_sub=esub[j]/hbarC;        // not s_an!!  convert to 1/fm^4 and do not rescale by factor!
        pj.u.x[0]=0;
        pj.u.x[1]=0;
        pj.eta_sigma = 1;
        pj.sigmaweight=stepx*stepy;
		pj.rho_weight = stepx*stepy;
        pj.Bulk = 0;
        pj.B=rBsub[j]*stepx*stepy;			// confirm with Jaki
        pj.S=rSsub[j]*stepx*stepy;			// confirm with Jaki
        pj.Q=rQsub[j]*stepx*stepy;			// confirm with Jaki
        pj.rhoB_an=rBsub[j];					// confirm with Jaki
        pj.rhoS_an=rSsub[j];					// confirm with Jaki
        pj.rhoQ_an=rQsub[j];					// confirm with Jaki
		pj.transverse_area = stepx*stepy;

		if (j==0)
		cout << "readICs_iccing(" << __LINE__ << "): "
			<< "SPH particles: "
			<< j << "   "
			<< pj.r.x[0] << "   " << pj.r.x[1] << "   "
			<< pj.e_sub << "   " << pj.rhoB_an << "   "
			<< pj.rhoS_an << "   " << pj.rhoQ_an << "   "
			<< pj.sigmaweight << endl;

		// make educated initial guess here for this particle's (T, mu_i) coordinates
		// (improve this in the future)
		pj.SPH_cell.T   = 500.0/197.3;	// rootfinder seems to work better going downhill than "uphill"
		pj.SPH_cell.muB = 0.0/197.3;
		pj.SPH_cell.muS = 0.0/197.3;
		pj.SPH_cell.muQ = 0.0/197.3;

        if (pj.e_sub>efcheck)	// impose freeze-out check for e, not s
        {
            pj.Freeze=0;
        }
        else
        {
            pj.Freeze=4;
            --kk;
            ++numpart;
        }
    }

    cout << "After freezeout (redundant): size = " << _Ntable3-numpart << endl;



}




}