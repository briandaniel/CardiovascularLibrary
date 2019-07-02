/*
 * SimpleMuscleFiber.cpp
 *
 *  Created on: Dec 28, 2018
 *      Author: brian
 */


#include "SimpleMuscleFiber.hpp"

double gaussTensionFunc( double lambda, double Ls0, double Lsmax, double Lsw)
{

	double LS;
	double geff;

	LS = Ls0*lambda;

	geff = exp( - pow( LS - Lsmax, 2) / (2*pow(Lsw,2) ) );

	return geff;

}



