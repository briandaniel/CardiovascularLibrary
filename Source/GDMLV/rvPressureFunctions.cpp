/*
 * rvPressureFunctions.cpp
 *
 *  Created on: Oct 21, 2017
 *      Author: brian
 */

#include "rvPressureFunctions.hpp"

void phiEllipse( double * coef, double nu0, double &phiMin, double &phiMax, double &phiMid )
{
	double c1, c2, a, b;
	double checkVal, rootVal;

	// pull values
	c1 = coef[0];
	c2 = coef[1];
	a = coef[2];
	b = coef[3];

	checkVal = pow( b, 2)*( 1 - pow(nu0 - c1 ,2) / pow( a, 2) );

	phiMid = c2;

	// if outside the ellipse, set impossible bounds for phi
	if ( checkVal < 0 )
	{

		phiMin = 2;
		phiMax = 1;

	}
	else
	{
		rootVal = sqrt( checkVal );
		phiMin = c2 - rootVal;
		phiMax = c2 + rootVal;

		// cout << "    " <<  c2 << "   " << phiMin << "   "  << phiMax << endl;
	}




}


bool checkPhiValue( double * coef, double nu0, double phi0 )
{

	double phiMin, phiMax, phiMid, phiShift, phi0Shifted;

	phiEllipse( coef, nu0, phiMin, phiMax, phiMid );

	phi0Shifted = mod(phi0, 2*PI);
	if ( phi0Shifted > (phiMid + PI) )
	{
		phi0Shifted = phi0Shifted- 2*PI;
	}
	else if (phi0Shifted < (phiMid - PI ) )
	{
		phi0Shifted = phi0Shifted + 2*PI;
	}

	bool inRegion = 0;

	if( phi0Shifted >= phiMin && phi0Shifted <= phiMax)
	{
		inRegion = 1;
	}


	// cout << "phiLow = " << phiMin << "  phiHigh = " << phiMax << "    phiMid = " << phiMid << endl;
	// cout << "   phi0 = " << phi0 << "   phi0Shifted = " << phi0Shifted <<  "   in region? " << inRegion << endl;



	return inRegion;



}


bool checkRVSurf( double * rvBoundaryCoef, double umu0, double nu0, double phi0 )
{
	bool rvSurfValue = 0;

	// must be on outer wall
	if (umu0 == 1)
	{
		rvSurfValue = checkPhiValue( rvBoundaryCoef, nu0, phi0 );
		//rvSurfValue = 1;
	}

	// cout << " nu0 = " << nu0 << "   phi0 = " << phi0 <<  "   umu0 " << umu0 << "   rvSurfValue " << rvSurfValue << "   phiCheckBool" << phiCheckBool << endl;

	return rvSurfValue;

}












