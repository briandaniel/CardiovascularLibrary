/*
 * LVDeformation.hpp
 *
 *  Created on: Sep 12, 2017
 *      Author: brian
 */

#ifndef LVDEFORMATION_HPP_
#define LVDEFORMATION_HPP_


// Standard constant definitions
#define PI 3.141592653589793
#define RTHIRD 1.259921049894873165 // 2^(1/3)
#define THIRD  0.333333333333333333 // 1/3

#include <math.h>
#include <complex>
#include <iostream>

void computeSevenParamDeformation( double muin0, double nu0, double phi0, double * q,
		double nuUp0, double fc0, double acube, double chmu0, double shmu0, double chmu0sq, double cnu0sq, double cphi0, double sphi0,
		double snu0, double snu0sq, double cnu0, double dfc0dnu0, double dfc0dphi0, double dmuin0dnu0, double dmuin0dphi0, double chmue,
		double &mu, double &nu, double &phi, double &shmu, double &chmu, double &shmusq, double &chmusq, double &snu, double &cnu, double &snusq, double &cnusq,
		double &fc, double &dfcdnu0, double &dfcdphi0, double &K, double &sinRatio, double & dnudnu0, double &dnudphi0, double &dphidnu0, double &dphidphi0,
		double &d2nudnu02, double &d2nudphi0dnu0, double & d2nudphi02, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0,
		double &dmudmu0, double &dmudnu0, double & dmudphi0);

void currentDeformationParameterSize( int & Nc, int & Nd, int & Ne, int & Nq );
void pullq_parameters( double *q, int Nc, int Nd, int Ne,  double *c, double *d, double *e);

double nuFunc( double nu0, double phi0, double nuUp0, double * e );
void nuDerivatives( double nu0, double phi0, double nuUp0, double * e, double &dnudnu0, double &dnudphi0 );
void nuSecondDerivatives( double nu0, double phi0, double nuUp0, double * e, double &d2nudnu02, double &d2nudphi0dnu0, double &d2nudphi02 );
double computeSinRatio( double nu0, double nu );

double phiFunc(  double nu0, double phi0, double nuUp0, double * d );
void phiDerivatives( double nu0, double phi0, double nuUp0, double * d, double &dphidnu0, double &dphidphi0 );
void phiSecondDerivatives( double nu0, double phi0, double nuUp0, double * d, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0 );

double ffunc( double * c, double acube, double cphi0, double sphi0, double nu0,
						  double fc0, double mue, double K, double sinRatio, double nu);
void fDerivatives( double acube, double sphi0, double cphi0, double nu0, double fc0, double dfc0dnu0, double dfc0dphi0, double R, double Rphi0,
		double Rnu0, double chmue, double snu, double cnu, double dnudphi0, double dnudnu0, double * c, double &dfcdnu0, double &dfcdphi0 );

double muFunc( double K, double sinRatio, double nu, double acube, double chmu0, double chmu0sq, double cnu0sq, double fc  );


using namespace std;

#endif /* LVDEFORMATION_HPP_ */
