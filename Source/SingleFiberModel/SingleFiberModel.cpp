/*
 * SingleFiberModel.cpp
 *
 *  Created on: Dec 11, 2018
 *      Author: brian
 */

#include "SingleFiberModel.hpp"



// Inputs: qin, qout, V, At
// Parameters: V0, Vw, Ta0, Tp0, cp, ls0, lsmax, lsw
// Outputs: P, dVdt, sigmaf
void singleFiberModelGaussTension( double qin, double qout, double V, double At, double V0, double Vw,  double Ta0, double Tp0,
		double cp, double cv, double v0, double ls0, double lsmax, double lsw, double & P, double & dVdt, double & sigmaf )
{

	double f, lambda, vs, h, sigmaa, sigmap;

	 // Volume conservation
	dVdt = qin - qout;

	// Fiber stretch in spherical model
	lambda = pow(  ( 1.0 + (3.0*V/Vw) )/( 1.0 + (3.0*V0/Vw) ), (1.0/3.0) );

	// Fiber shortening/lengthening velocity
	vs = -1.0/Vw*pow( 1.0 + (3.0*V)/Vw , -1 )*dVdt;

	// Active stress length dependence
	f = gaussTensionFunc( lambda, ls0, lsmax, lsw);

	// Velocity dependence
	h = hfuncFiber( vs, cv, v0 );

	// Active stress
	sigmaa = Ta0*At*f*h;

	// Passive stress
	sigmap = 0;
	if( lambda >= 1 )
	{
		sigmap = Tp0*( exp(cp*(lambda-1.0) ) - 1.0 );
	}

	// Total stress
	sigmaf = sigmaa + sigmap;

	// Pressure
	P = sigmaf/(1.0 + 3.0*V/Vw );


}


// Inputs: qin, qout, V, At
// Parameters: V0, Vw, Ta0, Tp0, cp, l0, la0, lam, lae, laf
// Outputs: P, dVdt, sigmaf
void singleFiberModel( double qin, double qout, double V, double At,
		double V0, double Vw, double Ta0, double Tp0, double cp, double cv,
		double v0, double l0, double la0, double lam, double lae, double laf,
		double & P, double & dVdt, double & sigmaf )
{

	double f, lambda, l, vs, h, sigmaa, sigmap;

	 // Volume conservation
	dVdt = qin - qout;

	// Fiber stretch in spherical model
	lambda = pow(  ( 1.0 + (3.0*V/Vw) )/( 1.0 + (3.0*V0/Vw) ), (1.0/3.0) );

	// Fiber length
	l = l0*lambda;

	// Fiber shortening/lengthening velocity
	vs = -1.0/Vw*pow( 1.0 + (3.0*V)/Vw , -1 )*dVdt;

	// Active stress length dependence
	f = ffuncFiber( l, l0, la0, lam, lae, laf );

	// Velocity dependence
	h = hfuncFiber( vs, cv, v0 );

	// Active stress
	sigmaa = Ta0*At*f*h;

	// Passive stress
	sigmap = 0;
	if( lambda >= 1 )
	{
		sigmap = Tp0*( exp(cp*(lambda-1.0) ) - 1.0 );
	}

	// Total stress
	sigmaf = sigmaa + sigmap;

	// Pressure
	P = sigmaf/(1.0 + 3.0*V/Vw );


}



double ffuncFiber( double l, double l0, double la0, double lam, double lae, double laf )
{
	double f = 0;

	if( la0 < l && l <= lam ){
		f = (l -  la0)/( lam -  la0);
	}else if(  lam < l && l <=  lae ){
		f = 1;
	}else if( l >  lae && l < laf ){
		f = (  laf - l )/(  laf -  lae ) ;
	}

	return f;
}

double hfuncFiber( double vs, double cv, double v0 )
{

	double h = ( 1 - (vs/v0) )/( 1 + cv*(vs/v0) );
	return h;

}







