/*
 * SingleFiberModel.hpp
 *
 *  Created on: Dec 11, 2018
 *      Author: brian
 */

#ifndef SINGLEFIBERMODEL_SINGLEFIBERMODEL_HPP_
#define SINGLEFIBERMODEL_SINGLEFIBERMODEL_HPP_

#include <math.h>
#include <iostream>
#include "../MuscleFiberFunctions/SimpleMuscleFiber.hpp"

void singleFiberModelGaussTension( double qin, double qout, double V, double At, double V0, double Vw,  double Ta0, double Tp0,
		double cp, double cv, double v0, double ls0, double lsmax, double lsw, double & P, double & dVdt, double & sigmaf );

void singleFiberModel( double qin, double qout, double V, double At,
		double V0, double Vw, double Ta0, double Tp0, double cp, double cv,
		double v0, double l0, double la0, double lam, double lae, double laf,
		double & P, double & dVdt, double & sigmaf );

double ffuncFiber( double l, double l0, double la0, double lam, double lae, double laf );

double hfuncFiber( double vs, double cv, double v0 );


#endif /* SINGLEFIBERMODEL_SINGLEFIBERMODEL_HPP_ */
