#ifndef TRANSPORT_COEFFICIENTS_H
#define TRANSPORT_COEFFICIENTS_H
#include <string>
#include <vector>
#include "settings.h"
#include "kernel.h"
#include "particle.h"
#include "new_matrix.h"

using std::string;
using std::vector;

class TransportCoeficients
{

  public:
    TransportCoeficients(){};
    ~TransportCoeficients(){};

    double getEta();
    double getZeta();
    // Matrix getKappa();
    double getTauShear();
    double getTauBulk();
    // Matrix getTauDiffusive();
    void setThermal(double currentT, double currentMuB,
                    double currentMuS, double currentMuQ);


  private:
    double T, muB, muS, muQ;
    double eta, zeta;
    // Matrix kappa;

    class constantEta
    {
      //move body to .cpp?
      constantEta(double inEta) : etaOVs(inEta) {}

    }



    









}



#endif
