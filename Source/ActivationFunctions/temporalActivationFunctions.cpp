/*
 * temporalActivationFunctions.cpp
 *
 *  Created on: Sep 12, 2018
 *      Author: brian
 */


#include "temporalActivationFunctions.hpp"




/* double simpleSinActivation( double t, double Tc, double Ta )
 *
 * Computes active curve during the last part of the cycle as a simple sine lobe
 * t: current time
 * Tc: cycle length
 * Ta: active duration
 *
 */
double simpleSinActivation( double t, double Tc, double Ta )
{

	double At = simpleSinActivationShifted( t, Tc, Ta, Tc-Ta );

	return At;

}


/* double simpleSinActivationShifted( double t, double Tc, double Ta, double Ts )
 *
 * Computes active curve during the last part of the cycle as a simple sine lobe
 * t: current time
 * Tc: cycle length
 * Ta: active duration
 * Ts: determines the shift of the activation away from the first part of the cycle
 * Note: it is OK to use negative numbers for Ts or go beyond the end of the cycle
 *
 */
double simpleSinActivationShifted( double t, double Tc, double Ta, double Ts )
{
	double At = 0;
	double t_cycle = mod(t-Ts,Tc);

    At = 0;
    if( t_cycle <= Ta )
    {
		At = sin(PI*t_cycle /Ta );
    }


	return At;

}



/* productActivation( double t, double Tc, double Ta, double Ts, double kt )
 *
 * Computes active curve during the last part of the cycle as a product
 * t: current time
 * Tc: cycle length
 * Ta: active duration
 * Ts: determines the shift of the activation away from the first part of the cycle
 * k: determines how wide the activation function is
 * Note: it is OK to use negative numbers for Ts or go beyond the end of the cycle
 *
 */
double productActivation( double t, double Tc, double Ta, double Ts, double kt )
{

	double At = 0;
	double t_cycle = mod(t-Ts,Tc);

    At = 0;
    if( t_cycle <= Ta )
    {
    	double strength = (tanh(kt) + 1)/2 + 0.05;
		double At1 = pow( sin(PI*t_cycle /Ta ), 4 );
		double At2 = exp( pow(t_cycle-Ta/2,2)/(strength*pow(Ta,2)));
		At = At1*At2;
    }


	return At;

}

/*
 * Computes active curve during the first part of the cycle but shifted by Ts as described in
 *
 * Mynard, J. P., et al. "A simple, versatile valve model for use in lumped parameter and one‐dimensional cardiovascular models."
 * International Journal for Numerical Methods in Biomedical Engineering 28.6-7 (2012): 626-641.
 *
 * t: current time
 * Tc: cycle length
 * Ta: active duration
 * Ts: determines the shift of the activation away from the first part of the cycle
 * m1: Contraction rate constant
 * m2: relaxation rate constant
 * tau1: systolic time constant
 * tau2: diastolic time constant
 *
 * hillMaxVal: PRECOMPUTED maximum value of the basic hill function which is scaled out
 *
 * Note: hillMaxVal MUST be PRECOMPUTED for this function to be normalized to 1
 *
 */
double twoHillActivation( double t, double m1, double m2, double tau1, double tau2, double Tc, double Ts, double hillMaxVal )
{
	double g1, g2, At;
	double t_cycle = mod(t-Ts,Tc);

    g1 = pow( t_cycle / tau1, m1 );
    g2 = pow( t_cycle / tau2, m2 );
    At = ( g1/(1 + g1)*(1/(1+g2)) )/hillMaxVal;

    return At;

}


// Auxiliary function for computing hillMaxVal
double twoHillActivationDerivative( double t, double m1, double m2, double tau1, double tau2, double Tc, double Ts )
{

	double t_cycle, g1, g2, g1prime, g2prime, AtPrime;

    t_cycle = mod(t - Ts, Tc);

    g1 = pow( t_cycle / tau1, m1 );
    g2 = pow( t_cycle / tau2, m2 );

	g1prime = (m1/tau1)*pow( t_cycle / tau1, m1-1);
	g2prime = (m2/tau2)*pow( t_cycle / tau2, m2-1);

	AtPrime = g1prime/( pow(1+g1,2)*(1+g2) ) - (g1*g2prime)/( (1+g1)*pow(1+g2,2) );

    return AtPrime;

}

// Computes the maximum of the un-normalized two hill activation function
double computeHillMax( double m1, double m2, double tau1, double tau2, double Tc )
{

	double Ts, tMax, hillMax, ta, tb, tc, t, hillk, hillaPrime, hillbPrime, hillcPrime;
	int k;
	double Nguess = 50;

    Ts = 0;
    tMax = 0;
    hillMax = -1e10;

    for ( k = 0; k < Nguess; k++)
    {
        t = Tc*(k/Nguess);
        hillk = twoHillActivation( t, m1, m2, tau1, tau2, Tc, Ts, 1.0 );
        if(hillk > hillMax)
        {
            hillMax = hillk;
            tMax = k/Nguess*Tc;
        }
    }


    ta = tMax - 1.0/Nguess*Tc;
    tb = tMax + 1.0/Nguess*Tc;

    hillaPrime = twoHillActivationDerivative( ta, m1, m2, tau1, tau2, Tc, Ts );
    hillbPrime = twoHillActivationDerivative( tb, m1, m2, tau1, tau2, Tc, Ts );

    // find critical point using bisection
    k = 0;
    while( max( fabs(hillaPrime), fabs(hillbPrime) ) > 1e-10 && k < 100)
    {

        tc = (ta+tb)/2.0;
        hillcPrime = twoHillActivationDerivative( tc, m1, m2, tau1, tau2, Tc, Ts );


        if( hillaPrime >= 0 && hillcPrime <= 0 )
        {
            tb = tc;
        }
        else
        {
            ta = tc;
        }

        hillaPrime = twoHillActivationDerivative( ta, m1, m2, tau1, tau2, Tc, Ts );
        hillbPrime = twoHillActivationDerivative( tb, m1, m2, tau1, tau2, Tc, Ts );

        k = k+1;
    }

    tMax = (ta+tb)/2.0;

    hillMax = twoHillActivation( tMax, m1, m2, tau1, tau2, Tc, Ts, 1.0 );

    return hillMax;
}

