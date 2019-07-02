/*
 * varyingElastance.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: brian
 */



#include <math.h>



/* 1. input variables: At: activation level, V: chamber volume,
 *    Pup: upstream chamber pressure: Pdown: downstream chamber pressure
 *
 * 2. input parameters: Emax: max elastance, kr: chamber flow resistance multiplier, Rup: upstream chamber flow resistance
 *    Vu: chamber unstressed volume,  ke: chamber exponential stiffness, P0: linear stiffness coef
 *
 * 3. output values: Pmax: max pressure generated, P: chamber pressure, Fi: flow in, Fo: flow out, dVdt: volume rate of change
 */

void computeVaryingElastanceModelUrsino( double At, double V, double Pup, double Pdown,
		double Emax, double kr, double Rup, double Vu, double ke, double P0,
		double & Pmax, double & P, double & Fi, double & Fo, double & dVdt);



/* Simpler varying elastance model taken from
 *
 * Mynard, J. P., et al. "A simple, versatile valve model for use in lumped parameter and one‐dimensional cardiovascular models."
 * International Journal for Numerical Methods in Biomedical Engineering 28.6-7 (2012): 626-641.
 *
 * 1. input variables: qin: flow into chamber, qout: flow out of chamber, V: chamber volume, At: activation level
 *
 * 2. input parameters: Emax: max elastance, Emin: min elastance, Ks: flow resistance coefficient, V0: zero pressure volume
 *
 * 3. output values: P: pressure, dVdt: volume rate of change
 */
void computeVaryingElastanceModelMynard( double qin, double qout, double V, double At,
		double Emax, double Emin, double Ks, double V0,
		double & P, double & dVdt );



/* computeVaryingElastanceModelAtrium0DCoupled(...)
 *
 * computes the varying elastance model from
 * Mynard, J. P., et al. "A simple, versatile valve model for use in lumped parameter and one‐dimensional cardiovascular models."
 * International Journal for Numerical Methods in Biomedical Engineering 28.6-7 (2012): 626-641.
 * but assuming that the inflow comes from a lumped parameter model with the specified properties, and solves for the pressure directly
 *
 * 1. input variables: qin: flow into chamber, qout: flow out of chamber, V: chamber volume, At: activation level
 *
 * 2. input parameters: Emax: max elastance, Emin: min elastance, Ks: flow resistance coefficient, V0: zero pressure volume
 *
 * 3. output values: P: pressure, dVdt: volume rate of change
 */
void computeVaryingElastanceModelAtrium0DCoupled( double Pup, double Rup, double qout, double V, double At,
		double Emax, double Emin, double Ks, double V0,
		double & P, double & dVdt );
