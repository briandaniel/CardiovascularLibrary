/*
 * ValveFunctions.cpp
 *
 *  Created on: Sep 13, 2018
 *      Author: brian
 */

#include "ValveFunctions.hpp"

/* void computeMynardValve(...)
 *
 * Computes the mynard valve formulated in
 *
 * Mynard, J. P., et al. "A simple, versatile valve model for use in lumped parameter and one‐dimensional cardiovascular models."
 * International Journal for Numerical Methods in Biomedical Engineering 28.6-7 (2012): 626-641.
 *
 * The inputs are: q: flow through the valve, zeta: the valve state, upP: the upstream pressure, downP: the downstream pressure
 * The parameters are:
 * rho: blood density, Aeff_max: max effective area of valve, Aeff_min: minimum effective valve area,
 * leff: effective distance blood travels through the valve, Kvo: valve opening rate, Kvc: Valve closing rate, d
 * eltaP_open: minimum pressure gradient for valve to open, deltaP_close: minimum pressure gradient for valve to close
 * The outputs are:
 * dqdt: flow rate of change, dZeta_dt: valve state rate of change
 *
 */
void computeMynardValve( double q, double zeta, double upP, double downP,
		double rho, double Aeff_max, double Aeff_min, double leff, double Kvo, double Kvc, double deltaP_open, double deltaP_close,
		double & dqdt, double & dZeta_dt )
{
	// local definitions
	double deltaP, Aeff, B, L;

    // Pressure difference is computed by subtracting the downstream pressure from the upstream pressure
    deltaP = upP - downP;

    // Compute the effective area
    Aeff = computeMynardEffectiveArea( zeta, Aeff_max, Aeff_min );

    // Flow determining parameters
    B = rho/(2 * pow(Aeff,2) );
    L = rho*leff/(Aeff);

    // Flow derivative
    dqdt = 1/L*( deltaP - B*q*fabs(q) );

    // Zeta derivative
    dZeta_dt = 0;
    if (deltaP > deltaP_open)
    {
        dZeta_dt = (1-zeta)*Kvo*(deltaP -  deltaP_open);
    }
    else if ( deltaP < deltaP_close )
    {
        dZeta_dt = zeta*Kvc*(deltaP - deltaP_close);
    }

}


// Computes the effective area based on the state variable
double computeMynardEffectiveArea( double zeta, double Aeff_max, double Aeff_min )
{
	double Aeff = zeta*( Aeff_max - Aeff_min) + Aeff_min;

	return Aeff;
}



// Computes the variable resistance valves
void valveResistances( double Pla, double Plv, double Ppao, double beta, double Rmvc, double Rmvo,
		double Raovc, double Raovo, double &Rmv, double &Raov )
{

	double denom1 = ( 1 + exp(- beta*(Pla - Plv) ) );
	double denom2 = ( 1 + exp(- beta*(Plv - Ppao) ) );

	// numerical safety
	// this is automatically computed by C++ anyways but its here to indicate the numerical issue that must be handled
	double mult1 = 1/denom1;
	double mult2 = 1/denom2;

	if( denom1 > 1e300)
	{
		mult1 = 0.0;
	}

	if( denom2 > 1e300 )
	{
		mult2 = 0.0;
	}

	Rmv = Rmvc -( Rmvc - Rmvo )*mult1;
	Raov = Raovc -( Raovc - Raovo )*mult2;

}


// Computes the derivatives of the variable resistance valves
void valveResistanceDerivatives( double Pla, double Plv, double Ppao, double beta, double Rmvc, double Rmvo,
		double Raovc, double Raovo,double &dRmv_dPlv, double &dRaov_dPlv, double &dRaov_dPpao )
{

	double denom1 = pow( exp( 0.5*beta*(Pla - Plv) ) +  exp( -0.5*beta*(Pla - Plv) ), 2);

	// numerical safety
	// this is automatically computed by C++ anyways but its here to indicate the numerical issue that must be handled
	double mult1 = beta/denom1;
	if ( denom1 > 1e300 )
		mult1 = 0;

	dRmv_dPlv = ( Rmvc - Rmvo )* mult1;

	double denom2 = pow( exp( 0.5*beta*(Plv - Ppao) ) +  exp( -0.5*beta*(Plv - Ppao) ), 2);

	// numerical safety
	// this is automatically computed by C++ anyways but its here to indicate the numerical issue that must be handled
	double mult2 = beta/denom2;
	if ( denom2 > 1e300 )
		mult2 = 0;

	dRaov_dPlv = -( Raovc - Raovo )* mult2;
	dRaov_dPpao = ( Raovc - Raovo )* mult2;

}








