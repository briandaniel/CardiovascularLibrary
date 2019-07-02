/*
 * cdmlv_node.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */



#include "cdmlv_node.hpp"

// This file contains the majority of the initial node setup


void CDMNode::setupNodeFromPrescribedLocation( double mu0In, double nu0In, double phi0In, ParamCDMLV & prm ){

	setPrescribedLocation( mu0In, nu0In, phi0In, prm );
	initialComputations( prm );

}


void CDMNode::setPrescribedLocation( double mu0In, double nu0In, double phi0In, ParamCDMLV & prm ){

	// these indices are not used if a prescribed location is given
	imu = -1;
	jnu = -1;
	kphi = -1;
	globalID = -1;
	localID = -1;

	// location
	mu0 = mu0In;
	nu0 = nu0In;
	phi0 = phi0In;

	// compute nuUp0
	nuUp0 = computeNuUp0( phi0, prm.nuUp0Spline );

	// Compute muin0, muout0
	muin0 = prm.endo0.evaluateSplines(nu0,phi0);
	muout0 = prm.epi0.evaluateSplines(nu0,phi0);

	// The umu0 coordinate ranges from [0,1] between muin0 and muout0
	umu0 = ( muout0 - mu0In )/( muout0 - muin0 );

	// Set whether or not this is on the endocardial surface if it is close to the boundary
	endoSurfNode = 0;
	if ( fabs(umu0) < 1e-10 ){
		endoSurfNode = 1;
	}

	// Set the local copies of parameter variables as needed
	a = prm.a;

}




void CDMNode::setupNodeFromLoadBalance(int globalID_in, int localID_in, ParamCDMLV & prm ){

	setIndicesFromLoadBalance(globalID_in, localID_in, prm );
	setLocationFromIndices( prm );
	initialComputations( prm );

}


// Set the indices from the load balance
void CDMNode::setIndicesFromLoadBalance(int globalID_in, int localID_in, ParamCDMLV & prm ){

	globalID = globalID_in;
	localID = localID_in;

	// location index
	imu = (int) floor( ( (double) globalID )/( prm.Nnu*prm.Nphi) );
	jnu = (int) floor( ( (double) ( globalID - imu*prm.Nnu*prm.Nphi))/( prm.Nphi) );
	kphi = globalID - imu*prm.Nnu*prm.Nphi - jnu*prm.Nphi;

}


// Set the location assuming that the indices have been set
void CDMNode::setLocationFromIndices( ParamCDMLV & prm ){

	// Compute phi0
    double dphi0 = 2 * PI / prm.Nphi;
	phi0 = dphi0*kphi;

	// Compute nu0
	nuUp0 = computeNuUp0( phi0, prm.nuUp0Spline );
	double dnu0 = (PI - nuUp0 ) / (prm.Nnu- 1);
	nu0 = nuUp0 + dnu0*jnu;

	// Compute muin0, muout0
	muin0 = prm.endo0.evaluateSplines(nu0,phi0);
	muout0 = prm.epi0.evaluateSplines(nu0,phi0);

	// The umu0 coordinate ranges from [0,1] between muin0 and muout0
	double dumu0 = ( 1.0 )/(prm.Nmu - 1.0);
	umu0 = imu*dumu0;

	// compute mu0
	mu0 = computeMu0( muin0, muout0, umu0 );

	// Set whether or not this is on the endocardial surface
	endoSurfNode = 0;
	if ( imu == 0 ){
		endoSurfNode = 1;
	}

	// Set the local copies of parameter variables as needed
	a = prm.a;

}



void CDMNode::initialComputations( ParamCDMLV & prm )
{
	// set to zero in case they are used without being needed
	vCrossTopNorm = 0.0;
	vCrossSideNorm = 0.0;

	acube = pow(a,3);

	prolate2xyz( mu0, nu0, phi0, a, x0, y0, z0);

	snu0 = sin(nu0);
	cnu0 = cos(nu0);

	snu0sq = pow(snu0,2);
	cnu0sq = pow(cnu0,2);

	shmu0 = sinh(mu0);
	chmu0 = cosh(mu0);

	shmu0sq = pow(shmu0,2);
	chmu0sq = pow(chmu0,2);

	gmu0 = a*sqrt( shmu0sq + snu0sq );
	gnu0 = gmu0;
	gphi0 = a*shmu0*snu0;

	gmult0 = gmu0*gnu0*gphi0;

	chmuin0 = cosh(muin0);
	chmuin0sq = pow(chmuin0,2);
	shmuin0 = sinh(muin0);

	cphi0 = cos(phi0);
	sphi0 = sin(phi0);

	// evaluate spline surface derivatives used in additional computations (derivative of the constant part is 0)
	prm.endo0.evaluateSplineDerivatives( nu0, phi0, dmuin0dnu0, dmuin0dphi0);
	prm.epi0.evaluateSplineDerivatives( nu0, phi0, dmuout0dnu0, dmuout0dphi0);

	computeFiberRotation( prm );

	// Compute the integration weights
	computeIntegrationWeights( prm );

	// compute the fc0 function and its derivatives
	fc0func( acube, chmuin0, chmuin0sq, shmuin0, dmuin0dnu0, dmuin0dphi0,
			 cnu0, cnu0sq, snu0, cphi0, fc0, dfc0dnu0, dfc0dphi0 );


}


void CDMNode::computeIntegrationWeights( ParamCDMLV & prm )
{
	double dmu0 = (muout0 - muin0) / (prm.Nmu - 1);

	double dnu0 = (PI - nuUp0 ) / (prm.Nnu- 1);

	double dphi0 = 2 * PI / prm.Nphi;


	double cimu0 = computeIntegralCoef( imu, prm.Nmu );
	double cjnu0 = computeIntegralCoef( jnu, prm.Nnu );
	double ckphi0 = computePeriodicIntegralCoef ( kphi, prm.Nphi );

	// The weight for a volume integral = (1.)(2.)(3.) with
	// 1. equal spacing deltas
	// 2. integration rule weights
	// 3. scale factors
	IV0weight = (dmu0*dnu0*dphi0)*(cimu0*cjnu0*ckphi0)*(gmu0*gnu0*gphi0);

	// The weight for a ( nu0 x phi0 ) area integral = (1.)(2.) with
	// 1. equal spacing deltas
	// 2. integration rule weights
	// note that there are no scale factors, because the scale is determined by the deformed configuration as well (computed later)
	IA0weight = (dnu0*dphi0)*(cjnu0*ckphi0);

	// Same thing for the top integral
	IA0weight_top = (dmu0*dphi0)*(cimu0*ckphi0);


}



