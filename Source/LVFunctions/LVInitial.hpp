/*
 * LVInitial.hpp
 *
 *  Created on: Sep 12, 2017
 *      Author: brian
 */

#ifndef LVINITIAL_HPP_
#define LVINITIAL_HPP_


// Standard constant definitions
#define PI 3.141592653589793
#define RTHIRD 1.259921049894873165 // 2^(1/3)
#define THIRD  0.333333333333333333 // 1/3


#include <math.h>
#include <complex>
#include "CubicSplines1D/cubicSplinesFunctions1D.hpp"

double computeNuUp0( double phi0, Spline1D & nuUp0Spline );
double computeMu0( double muin0, double muout0, double umu0 );

void fc0func( double acube, double chmuin0, double chmuin0sq, double shmuin0,  double dmuin0dnu0, double dmuin0dphi0,
		double cnu0, double cnu0sq, double snu0, double cphi0, double & fc0, double & dfc0dnu0, double & dfc0dphi0 );

using namespace std;

#endif /* LVINITIAL_HPP_ */
