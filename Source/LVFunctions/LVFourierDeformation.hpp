/*
 * LVFourierDeformation.hpp
 *
 *  Created on: Jun 17, 2018
 *      Author: brian
 */

#ifndef LVFUNCTIONS_LVFOURIERDEFORMATION_HPP_
#define LVFUNCTIONS_LVFOURIERDEFORMATION_HPP_



// Standard constant definitions
#define PI 3.141592653589793
#define RTHIRD 1.259921049894873165 // 2^(1/3)
#define THIRD  0.333333333333333333 // 1/3

// standard includes
#include <math.h>
#include <complex>
#include <iostream>

#include "UtilityFunctions/utilityFunctions.hpp"


// namespace
using namespace std;

// local classes

class FourierDeformation {

  public:

	// size of fourier deformations
	int muinNnuBasis, muinNphiBasis;
	int nuNnuBasis, nuNphiBasis;
	int phiNnuBasis, phiNphiBasis;

	int Nc, Nd, Ne;

	// main function that computes several values required by the LV model
	void computeDeformation( double muin0, double nu0, double phi0, double * q,
		double nuUp0, double fc0, double acube, double chmu0, double shmu0, double chmu0sq, double cnu0sq, double cphi0, double sphi0,
		double snu0, double snu0sq, double cnu0, double dfc0dnu0, double dfc0dphi0, double dmuin0dnu0, double dmuin0dphi0,
		double &mu, double &nu, double &phi, double &shmu, double &chmu, double &shmusq, double &chmusq, double &snu, double &cnu, double &snusq, double &cnusq,
		double &fc, double &dfcdnu0, double &dfcdphi0, double &K, double &sinRatio, double & dnudnu0, double &dnudphi0, double &dphidnu0, double &dphidphi0,
		double &d2nudnu02, double &d2nudphi0dnu0, double & d2nudphi02, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0,
		double &dmudmu0, double &dmudnu0, double & dmudphi0 );

	// local functions
	void setSize( int muinNnuBasisIn, int muinNphiBasisIn, int nuNnuBasisIn,
			int nuNphiBasisIn, int phiNnuBasisIn, int phiNphiBasisIn );
	void computeNParams( int & NqOut, int & NcOut, int & NdOut, int & NeOut );
	void pullq_parameters( double *q, double *coef_cf, double *coef_df, double *coef_ef );
	void printNModesToConsole();

	void nuFunc( double nu0, double phi0, double * ef, double nuUp0,
				double & nu, double & dnudnu0, double & dnudphi0, double & d2nudnu02, double & d2nudphi0dnu0, double & d2nudphi02 );

	void phiFunc( double nu0, double phi0, double * df, double nuUp0,
				double & phi, double & dphidnu0, double & dphidphi0, double & d2phidnu02, double & d2phidphi0dnu0, double & d2phidphi02 );

	void muinFunc( double muin0, double nu0, double phi0, double * cf, double dmuin0dnu0, double dmuin0dphi0, double nuUp0,
				double & muin, double & dmuindnu0, double & dmuindphi0 );

	// constructor/destructor
	FourierDeformation();
	FourierDeformation( int muinNnuBasisIn, int muinNphiBasisIn, int nuNnuBasisIn,
				int nuNphiBasisIn, int phiNnuBasisIn, int phiNphiBasisIn );

};

// local functions
/*
double computeFourierValue( double nu0, double phi0, double * coef, int Nnu0, int Nphi0 );
void computeFourierValueDerivatives( double nu0, double phi0, double * coef, int Nnu0, int Nphi0, double & dfdnu0, double & dfdphi0 );
void computeFourierValueSecondDerivatives( double nu0, double phi0, double * coef, int Nnu0, int Nphi0,
		double &d2fdnu02, double &d2fdphi0dnu0, double &d2fdphi02 );
*/
void computeFourierValueAndDerivatives( double nu0, double phi0, double * coef, int Nnu0, int Nphi0, double nuUp0,
		double & f, double & dfdnu0, double & dfdphi0, double &d2fdnu02, double &d2fdphi0dnu0, double &d2fdphi02 );

void computeFourierFunctionValuesAndDerivativesSinglePoint( double nu0, double phi0, int i, int jcos, int jsin, double nuUp0,
		double & fnu0, double &  dfnu0_dnu0, double & d2fnu0_dnu02,
		double & fcphi0, double & dfcphi0_dphi0, double & d2fcphi0_dphi02,
		double & fsphi0, double & dfsphi0_dphi0, double & d2fsphi0_dphi02 );

double computeSinRatioFourier( double nu0, double nu );

double fcfuncFourier( double muin, double fc0, double acube, double R, double cnusq);

double muFuncFourier( double K, double sinRatio, double nu, double acube, double chmu0, double chmu0sq, double cnu0sq, double fc  );

void fcDerivativesFourier( double  muin, double dmuindnu0, double dmuindphi0,
		double dfc0dnu0, double dfc0dphi0, double acube, double R, double Rnu0, double Rphi0, double cnu,
		double cnusq, double snu, double dnudnu0, double dnudphi0, double & dfcdnu0, double & dfcdphi0 );

#endif /* LVFUNCTIONS_LVFOURIERDEFORMATION_HPP_ */
