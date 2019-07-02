/*
 * temporalActivationFunctions.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: brian
 */

#ifndef ACTIVATIONFUNCTIONS_TEMPORALACTIVATIONFUNCTIONS_HPP_
#define ACTIVATIONFUNCTIONS_TEMPORALACTIVATIONFUNCTIONS_HPP_



#include <math.h>
#include "UtilityFunctions/utilityFunctions.hpp"

#define PI 3.141592653589793

// local functions
/* double simpleSinActivation( double t, double Tc, double Ta )
 *
 * Computes active curve during the last part of the cycle as a simple sine lobe
 * t: current time
 * Tc: cycle length
 * Ta: active duration
 *
 */
double simpleSinActivation( double t, double Tc, double Ta );


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
double simpleSinActivationShifted( double t, double Tc, double Ta, double Ts );



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
double productActivation( double t, double Tc, double Ta, double Ts, double k );

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
double twoHillActivation( double t, double m1, double m2, double tau1, double tau2, double Tc, double Ts, double hillMaxVal );

// Auxiliary function for computing hillMaxVal
double twoHillActivationDerivative( double t, double m1, double m2, double tau1, double tau2, double Tc, double Ts );

// Computes the maximum of the un-normalized two hill activation function
double computeHillMax( double m1, double m2, double tau1, double tau2, double Tc );



#endif /* ACTIVATIONFUNCTIONS_TEMPORALACTIVATIONFUNCTIONS_HPP_ */










