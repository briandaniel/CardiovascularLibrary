/*
 * ValveFunctions.hpp
 *
 *  Created on: Sep 13, 2018
 *      Author: brian
 */

#ifndef VALVEFUNCTIONS_HPP_
#define VALVEFUNCTIONS_HPP_

#include <math.h>

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
		double & dqdt, double & dZeta_dt );

// Computes the effective area based on the state variable
double computeMynardEffectiveArea( double zeta, double Aeff_max, double Aeff_min );


// Computes the variable resistance valves
void valveResistances( double Pla, double Plv, double Ppao, double beta, double Rmvc, double Rmvo,
		double Raovc, double Raovo, double &Rmv, double &Raov );

// Computes the derivatives of the variable resistance valves
void valveResistanceDerivatives( double Pla, double Plv, double Ppao, double beta, double Rmvc, double Rmvo,
		double Raovc, double Raovo,double &dRmv_dPlv, double &dRaov_dPlv, double &dRaov_dPpao );





#endif /* VALVEFUNCTIONS_HPP_ */
