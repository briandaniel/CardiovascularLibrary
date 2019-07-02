/*
 * deformation.cpp
 *
 *  Created on: Sep 12, 2017
 *      Author: brian
 */

#include "../LVFunctions/LVInitial.hpp"



double computeNuUp0( double phi0, Spline1D & nuUp0Spline )
{
	phi0 = mod(phi0, 2*PI);

	double nuUp0 = nuUp0Spline.evaluateSpline( phi0 );

	return nuUp0;
}



double computeMu0( double muin0, double muout0, double umu0 )
{
	double mu0;

	mu0 = muin0*(1.0 - umu0) + muout0*umu0;

	return mu0;
}


void fc0func( double acube, double chmuin0, double chmuin0sq, double shmuin0,  double dmuin0dnu0, double dmuin0dphi0,
		double cnu0, double cnu0sq, double snu0, double cphi0, double & fc0, double & dfc0dnu0, double & dfc0dphi0 )
{
	double third = 1.0/3.0;

	fc0 = acube*( chmuin0*( third*chmuin0sq - cnu0sq ) - ( third - cnu0sq ) ) ;
	dfc0dnu0 = acube*( chmuin0sq*shmuin0*dmuin0dnu0 - cnu0sq*shmuin0*dmuin0dnu0 + 2*cnu0*snu0*(chmuin0 - 1) );
	dfc0dphi0 = acube*( chmuin0sq*shmuin0*dmuin0dphi0 - cnu0sq*shmuin0*dmuin0dphi0 );

}


