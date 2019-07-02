/*
 * LVSplineDeformation.hpp
 *
 *  Created on: Jun 5, 2018
 *      Author: brian
 */

#ifndef LVFUNCTIONS_LVBICUBICDEFORMATION_HPP_
#define LVFUNCTIONS_LVBICUBICDEFORMATION_HPP_



// Standard constant definitions
#define PI 3.141592653589793
#define RTHIRD 1.259921049894873165 // 2^(1/3)
#define THIRD  0.333333333333333333 // 1/3

// standard includes
#include <math.h>
#include <complex>
#include <iostream>

// local includes
#include "ProlateSplines/prolateSplines.hpp"

// namespace
using namespace std;

// local classes

class BicubicDeformation {

  public:

	// main function that computes several values required by the LV model
	void computeDeformation( double muin0, double nu0, double phi0,
		double nuUp0, double fc0, double acube, double chmu0, double shmu0, double chmu0sq, double cnu0sq, double cphi0, double sphi0,
		double snu0, double snu0sq, double cnu0, double dfc0dnu0, double dfc0dphi0, double dmuin0dnu0, double dmuin0dphi0,
		double &mu, double &nu, double &phi, double &shmu, double &chmu, double &shmusq, double &chmusq, double &snu, double &cnu, double &snusq, double &cnusq,
		double &fc, double &dfcdnu0, double &dfcdphi0, double &K, double &sinRatio, double & dnudnu0, double &dnudphi0, double &dphidnu0, double &dphidphi0,
		double &d2nudnu02, double &d2nudphi0dnu0, double & d2nudphi02, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0,
		double &dmudmu0, double &dmudnu0, double & dmudphi0 );

	// local functions
	void setupSplines( int muinNnuGrid, int muinNphiGrid, int nuNnuGrid, int nuNphiGrid,
			int phiNnuGrid, int phiNphiGrid, double nuUp0Min, double muin0_const );
	void updateBicubicDeformationSplines( double * q, int Nc, int Nd, int Ne );


	// constructor/destructor
	BicubicDeformation();
	BicubicDeformation( int muinNnuGrid, int muinNphiGrid, int nuNnuGrid, int nuNphiGrid,
			int phiNnuGrid, int phiNphiGrid, double nuUp0Min, double muin0_const );
	~BicubicDeformation();


	ProlateSplines nuSpline;
	ProlateSplines phiSpline;
	ProlateSplines fmuinSpline;


};

// local functions
void pullq_parameters_bicub( double *q, int Nc, int Nd, int Ne,  double *sc, double *sd, double *se);

double nuFuncBicub( double nu0, double phi0, ProlateSplines & nuSpline);
void nuDerivativesBicub( double nu0, double phi0, ProlateSplines & nuSpline, double &dnudnu0, double &dnudphi0 );
void nuSecondDerivativesBicub( double nu0, double phi0, ProlateSplines & nuSpline, double &d2nudnu02, double &d2nudphi0dnu0, double &d2nudphi02 );
double computeSinRatioBicub( double nu0, double nu );

double phiFuncBicub(  double nu0, double phi0, ProlateSplines & phiSpline);
void phiDerivativesBicub( double nu0, double phi0, ProlateSplines & phiSpline, double &dphidnu0, double &dphidphi0 );
void phiSecondDerivativesBicub( double nu0, double phi0, ProlateSplines & phiSpline, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0 );

double fcfuncBicub( double muin0, double nu0, double phi0, ProlateSplines & fmuinSpline, double fc0, double acube, double R, double cnusq);
double muFuncBicub( double K, double sinRatio, double nu, double acube, double chmu0, double chmu0sq, double cnu0sq, double fc  );

void fcDerivativesBicub(double muin0, double nu0, double phi0, ProlateSplines & fmuinSpline, double dfc0dnu0, double dfc0dphi0, double acube,
		double dmuin0dnu0, double dmuin0dphi0, double R, double Rnu0, double Rphi0, double cnu, double cnusq, double snu, double dnudnu0, double dnudphi0,
		double & dfcdnu0, double & dfcdphi0 );

#endif /* LVFUNCTIONS_LVBICUBICDEFORMATION_HPP_ */
