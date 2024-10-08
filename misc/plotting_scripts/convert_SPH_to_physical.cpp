#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

constexpr double PI 	= 3.1415926535897932384626433;
constexpr double h 		= 0.3;
constexpr double knorm 	= 10.0/(7.*PI*h*h);

double kernel( double x, double y )
{
	double r = sqrt(x*x+y*y);
	double q = r/h;

	if(q>=2)
		return 0;
	else if(q>=1)
		return 0.25*knorm*(2.0-q)*(2.0-q)*(2.0-q);
	else
	{
		double qq = q*q;
		return knorm*(1.0 - 1.5*qq + 0.75*q*qq);
	}
}

int nX = 0, nY = 0;

int grid_index(int ix, int iy)
{
  return ix * nY + iy;
}


int main( int argc, char ** argv )
{
	string resultsDirectory 	= "./results";
	string fileStemToReadIn 	= argv[1];
	string fileStemToPrintTo 	= argv[2];
	string indexToProcess 		= argv[3];

	string infilename = resultsDirectory + "/" + fileStemToReadIn + indexToProcess + ".dat";
	string outfilename = resultsDirectory + "/" + fileStemToPrintTo + indexToProcess + ".dat";

	vector<double> xvec, yvec, pvec, Tvec, muBvec, muQvec, muSvec, evec, Bvec, Svec, Qvec, svec;
	double tau, x, y, pin, Tin, muBin, muSin, muQin, ein, Bin, Sin, Qin, sin;

	// read in data here
	int linecount = 0;
	long long nSPH = 0;
	double dummy = 0.0;
	ifstream infile( infilename.c_str() );
	if (infile.is_open())
	{
		string line;
		while ( getline (infile, line) )
		{
			istringstream iss(line);
			if ( linecount++ < 1 )						// skip first row
				continue;
			else
      {
				iss >> dummy >> tau >> x >> y >> pin
					>> Tin >> muBin >> muSin >> muQin
					>> ein >> Bin >> Sin >> Qin >> sin;	// discard rest of line after this
        while (iss) iss >> dummy;             // except for the last argument
      }

//cout << "CHECK: " << dummy << endl;
//exit(1);

      if (int(dummy) < 4 )
      {
        xvec.push_back( x );
        yvec.push_back( y );
        pvec.push_back( pin );
        Tvec.push_back( Tin );
        muBvec.push_back( muBin );
        muSvec.push_back( muSin );
        muQvec.push_back( muQin );
        evec.push_back( ein );
        Bvec.push_back( Bin );
        Svec.push_back( Sin );
        Qvec.push_back( Qin );
        svec.push_back( sin );
        nSPH++;
      }
		}
	}

	infile.close();

	// now loop through and compute physical quantities
  constexpr double TINY = 1e-10;
	const double dx = 0.05, dy = 0.05;
	const double xmin = -15.0, ymin = -15.0;
	const double xmax = -xmin, ymax = -ymin;

	/*
	ofstream outfile( outfilename.c_str() );

  for (double x_local = xmin; x_local < xmax + 1e-10; x_local += dx )
	for (double y_local = ymin; y_local < ymax + 1e-10; y_local += dy )
	{
		double normalization 				= 1e-100;	// protects from dividing by zero below
		double temperature 					= 0.0;
		double baryon_chemical_potential 	= 0.0;
		double strange_chemical_potential 	= 0.0;
		double electric_chemical_potential 	= 0.0;
		double energy_density 				= 0.0;
		double baryon_density 				= 0.0;
		double strange_density 				= 0.0;
		double electric_density 			= 0.0;
		double entropy_density 				= 0.0;

		// loop over SPH particles
		for (int iSPH = 0; iSPH < nSPH; iSPH++)
		{
      const double delta_x = x_local - xvec[iSPH];
      const double delta_y = y_local - yvec[iSPH];
      if (delta_x*delta_x+delta_y*delta_y > 4.0*h*h) continue;

			double kern 				 = kernel(delta_x, delta_y);
			normalization 				+= kern;
			energy_density 				+= kern * evec[iSPH];
			baryon_density 				+= kern * Bvec[iSPH];
			strange_density 			+= kern * Svec[iSPH];
			electric_density 			+= kern * Qvec[iSPH];
			temperature 				+= kern * Tvec[iSPH];
			baryon_chemical_potential 	+= kern * muBvec[iSPH];
			strange_chemical_potential 	+= kern * muSvec[iSPH];
			electric_chemical_potential += kern * muQvec[iSPH];
			entropy_density 			+= kern * svec[iSPH];
		}

		// normalize results
		energy_density 				/= normalization;
		baryon_density 				/= normalization;
		strange_density 			/= normalization;
		electric_density 			/= normalization;
		temperature 				/= normalization;
		baryon_chemical_potential 	/= normalization;
		strange_chemical_potential 	/= normalization;
		electric_chemical_potential /= normalization;
		entropy_density 			/= normalization;

		outfile << setw(12) << setprecision(8) << scientific
			<< tau << "   " << x_local << "   " << y_local << "   "
			<< temperature << "   " << baryon_chemical_potential << "   "
			<< strange_chemical_potential << "   " << electric_chemical_potential << "   "
			<< energy_density << "   " << baryon_density << "   " << strange_density << "   "
			<< electric_density << "   " << entropy_density << "\n";

	}*/

  vector<double> xGrid, yGrid;
  for (double x_local = xmin; x_local < xmax + 1e-10; x_local += dx )
    xGrid.push_back( x_local );
  for (double y_local = ymin; y_local < ymax + 1e-10; y_local += dy )
    yGrid.push_back( y_local );

  nX = xGrid.size();
  nY = yGrid.size();

  vector<double> normGrid(nX*nY, 1e-100);
  vector<double> TGrid(nX*nY);
  vector<double> muBGrid(nX*nY);
  vector<double> muSGrid(nX*nY);
  vector<double> muQGrid(nX*nY);
  vector<double> eGrid(nX*nY);
  vector<double> BGrid(nX*nY);
  vector<double> SGrid(nX*nY);
  vector<double> QGrid(nX*nY);
  vector<double> sGrid(nX*nY);

  // loop over SPH particles
  for (int iSPH = 0; iSPH < nSPH; iSPH++)
  {
    const double x0  = xvec[iSPH];
    const double y0  = yvec[iSPH];
    const double e   = evec[iSPH];
    const double B   = Bvec[iSPH];
    const double S   = Svec[iSPH];
    const double Q   = Qvec[iSPH];
    const double T   = Tvec[iSPH];
    const double muB = muBvec[iSPH];
    const double muS = muSvec[iSPH];
    const double muQ = muQvec[iSPH];
    const double s   = svec[iSPH];

    const int ix0 = int((x0-xmin)/dx);
    const int iy0 = int((y0-ymin)/dy);
    const int x_box_size = 3.0*h/dx;
    const int y_box_size = 3.0*h/dy;
    const int ix_min = std::max(0, ix0-x_box_size);
    const int ix_max = std::min((int)xGrid.size()-1, ix0+x_box_size);
    const int iy_min = std::max(0, iy0-y_box_size);
    const int iy_max = std::min((int)yGrid.size()-1, iy0+y_box_size);

    for (int ix = ix_min; ix <= ix_max; ix++)
    for (int iy = iy_min; iy <= iy_max; iy++)
    {
      int index        = grid_index(ix, iy);
      double delta_x   = xGrid[ix] - x0;
      double delta_y   = yGrid[iy] - y0;
      double kern      = kernel(delta_x, delta_y);

			normGrid[index] += kern;
			eGrid[index]    += kern * e;
			BGrid[index]    += kern * B;
			SGrid[index]    += kern * S;
			QGrid[index]    += kern * Q;
			TGrid[index]    += kern * T;
			muBGrid[index]  += kern * muB;
			muSGrid[index]  += kern * muS;
			muQGrid[index]  += kern * muQ;
			sGrid[index]    += kern * s;
    }
  }

	ofstream outfile( outfilename.c_str() );

  for (int ix = 0; ix < xGrid.size(); ix++)
  for (int iy = 0; iy < yGrid.size(); iy++)
  {
    int index       = grid_index(ix, iy);
    double x0       = xGrid[ix];
    double y0       = yGrid[iy];

    // normalize results
    double norm     = normGrid[index];
		eGrid[index]   /= norm;
		BGrid[index]   /= norm;
		SGrid[index]   /= norm;
		QGrid[index]   /= norm;
		TGrid[index]   /= norm;
		muBGrid[index] /= norm;
		muSGrid[index] /= norm;
		muQGrid[index] /= norm;
		sGrid[index]   /= norm;

		outfile << setw(12) << setprecision(8) << scientific
			<< tau << "   "   << x0 << "   " << y0 << "   "
			<< TGrid[index]   << "   " << muBGrid[index] << "   "
			<< muSGrid[index] << "   " << muQGrid[index] << "   "
			<< eGrid[index]   << "   " << BGrid[index]   << "   " << SGrid[index] << "   "
			<< QGrid[index]   << "   " << sGrid[index]   << "\n";
  }

	outfile.close();

	return 0;
}
