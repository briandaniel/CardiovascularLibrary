/*
 * rvPressureFunctions.hpp
 *
 *  Created on: Oct 21, 2017
 *      Author: brian
 */

#ifndef SOURCE_RVPRESSUREFUNCTIONS_HPP_
#define SOURCE_RVPRESSUREFUNCTIONS_HPP_


#include "UtilityFunctions/utilityFunctions.hpp"

void phiEllipse( double * coef, double nu0, double &phiMin, double &phiMax, double &phiMid );
bool checkPhiValue( double * coef, double nu0, double phi0 );
bool checkRVSurf( double * rvBoundaryCoef, double umu0, double nu0, double phi0 );

#endif /* SOURCE_RVPRESSUREFUNCTIONS_HPP_ */
