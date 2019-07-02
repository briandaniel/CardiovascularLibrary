/*
 * LVReferenceAdjustments.cpp
 *
 *  Created on: May 16, 2019
 *      Author: brian
 */


#include "LVReferenceAdjustments.hpp"

double LVReferenceShape::focalLength(){

	double a = a_original + delta_a;

	return a;

}

void LVReferenceShape::computeRefSurfaceValuesOnly( double nu0, double phi0, double & muin0, double & muout0 )
{
	double dmuin0dnu0, dmuin0dphi0, dmuout0dnu0,dmuout0dphi0;
	computeRefSurfaceValuesAndDerivatives(nu0, phi0,
			 muin0,  muout0,  dmuin0dnu0,  dmuin0dphi0,  dmuout0dnu0,  dmuout0dphi0);
}

void LVReferenceShape::computeRefSurfaceValuesAndDerivatives( double nu0, double phi0,
		double & muin0, double & muout0, double & dmuin0dnu0, double & dmuin0dphi0, double & dmuout0dnu0, double & dmuout0dphi0 )
{

	muin0 = endo0.evaluateSplines(nu0,phi0);
	muout0 = epi0.evaluateSplines(nu0,phi0);

	endo0.evaluateSplineDerivatives( nu0, phi0, dmuin0dnu0, dmuin0dphi0);
	epi0.evaluateSplineDerivatives( nu0, phi0, dmuout0dnu0, dmuout0dphi0);

	if( useRefAdj )
	{
		double muin0s, muout0s, dmuin0sdnu0, dmuin0sdphi0, dmuout0sdnu0, dmuout0sdphi0;
		computeRefSurfaceValuesAdjustment( nu0, phi0,
			 muin0, muout0, dmuin0dnu0, dmuin0dphi0, dmuout0dnu0, dmuout0dphi0,
			 muin0s, muout0s, dmuin0sdnu0, dmuin0sdphi0, dmuout0sdnu0, dmuout0sdphi0 );


		//cout << muin0 << " --> " <<  muin0s << endl;
		//print1DVector(q0s);
	    // Replace the normal values with the adjusted values
	    muin0 = muin0s;
	    muout0 = muout0s;
	    dmuin0dnu0 = dmuin0sdnu0;
	    dmuin0dphi0 = dmuin0sdphi0;
	    dmuout0dnu0 = dmuout0sdnu0;
	    dmuout0dphi0 = dmuout0sdphi0;

	}


}


void LVReferenceShape::computeRefSurfaceValuesAdjustment( double nu0, double phi0,
		double muin0, double muout0, double dmuin0dnu0, double dmuin0dphi0, double dmuout0dnu0, double  dmuout0dphi0,
		double & muin0s, double & muout0s, double & dmuin0sdnu0, double & dmuin0sdphi0, double & dmuout0sdnu0, double & dmuout0sdphi0 )
{
	//------------------------- 1. compute muin0 -----------------------------//

	// Compute the adjustment values at the endocardial surface
	double fs, dfsdnu0, dfsdphi0, d2fsdnu02, d2fsdphi0dnu0, d2fsdphi02;
	vector<double> q0s_fourier( q0s.size() - 1, 0);
	for(int i = 0; i < q0s_fourier.size(); i++)
	{
		q0s_fourier[i] = q0s[i+1];
	}
	computeFourierValueAndDerivatives( nu0, phi0, q0s_fourier.data(), refAdjust_muinNnuBasis, refAdjust_muinNphiBasis, nuUp0Min,
			fs, dfsdnu0, dfsdphi0, d2fsdnu02, d2fsdphi0dnu0, d2fsdphi02 );

	muin0s = muin0 + fs;
	dmuin0sdnu0 = dmuin0dnu0 + dfsdnu0;
	dmuin0sdphi0 = dmuin0dphi0 + dfsdphi0;


	//------------------------- 2. determine muout0 based on volume conservation -----------------------------//

	// compute the h values and derivatives (these are dummies for computational convenience)
	double h_muin0, dhdnu0_muin0, dhdphi0_muin0;
	double h_muout0, dhdnu0_muout0, dhdphi0_muout0;
	double h_muin0s, dhdnu0_muin0s, dhdphi0_muin0s;

	hFuncAndDerivatives_refFuncs( muin0, nu0, dmuin0dnu0, dmuin0dphi0, h_muin0, dhdnu0_muin0, dhdphi0_muin0);
	hFuncAndDerivatives_refFuncs( muout0, nu0, dmuout0dnu0, dmuout0dphi0, h_muout0, dhdnu0_muout0, dhdphi0_muout0);
	hFuncAndDerivatives_refFuncs( muin0s, nu0, dmuin0sdnu0, dmuin0sdphi0, h_muin0s, dhdnu0_muin0s, dhdphi0_muin0s);

	// compute another dummy (gs) and derivatives
	double gs, dgsdnu0, dgsdphi0;
	double as = a_original + delta_a;
	double afactor = pow(a_original/as,3);
	gs = afactor*( h_muout0 - h_muin0 )+ h_muin0s;
	dgsdnu0 = afactor*( dhdnu0_muout0 - dhdnu0_muin0 ) + dhdnu0_muin0s;
	dgsdphi0 = afactor*( dhdphi0_muout0 - dhdphi0_muin0 ) + dhdphi0_muin0s;

	// solve for mu
	double s1, s2, r3;
	std::complex <double> root3, coshMuComplex;
	s1 = - pow( cos(nu0), 2);
	s2 = - gs;
	root3 = pow( sqrt( (std::complex<double>)  ( 4*pow(s1,3) + 9*pow(s2,2) ) ) - 3*s2, 1.0/3.0 );
	coshMuComplex = root3/RTHIRD - RTHIRD*s1/root3;
	muout0s = acosh( coshMuComplex.real() );

	// compute muout0s derivatives
	dmuout0sdnu0 = ( dgsdnu0 - 2.0*cosh(muout0s)*cos(nu0)*sin(nu0) )/( sinh(muout0s)*( pow( cosh(muout0s),2 ) - pow( cos(nu0), 2 ) ) );
	dmuout0sdphi0 = ( dgsdphi0 )/( sinh(muout0s)*( pow( cosh(muout0s), 2 ) - pow( cos(nu0), 2)  ) );

	// in case of over contraction simply set mu = 0
	if(muout0s != muout0s)
	{
		muout0s = 0.0;
	}



}


void hFuncAndDerivatives_refFuncs( double mu0, double nu0, double dmu0dnu0, double dmu0dphi0,
		double & h, double & dhdnu0, double & dhdphi0 )
{
	double chmu0 = cosh(mu0);
	double chmu0sq = pow( cosh(mu0), 2 );
	double cnu0sq = pow( cos(nu0), 2 );

	h = chmu0*( THIRD*chmu0sq - cnu0sq );

	dhdnu0 = dmu0dnu0*sinh(mu0)*( chmu0sq - cnu0sq ) + 2*cos(nu0)*sin(nu0)*cosh(mu0);
	dhdphi0 = dmu0dphi0*sinh(mu0)*( chmu0sq - cnu0sq );

}




void LVReferenceShape::setRefSurfAdjValues( vector<double> & q0sIn )
{
	// The first value adjusts the ellipse axis
	delta_a = q0sIn[0];

	// the remaining values adjust walls
	q0s.resize(q0sIn.size());
	for (int i = 0; i < q0sIn.size(); i++ )
	{
		q0s[i] = q0sIn[i];
	}

}


void LVReferenceShape::setupRefSurfAdj(  int refAdjust_muinNnuBasisIn, int refAdjust_muinNphiBasisIn )
{
	// Set up the reference adjustment fourier deformation
	refAdjust_muinNnuBasis = refAdjust_muinNnuBasisIn;
	refAdjust_muinNphiBasis = refAdjust_muinNphiBasisIn;

	// there are no adjustments to nu or phi
	int nuNnuBasis = 0;
	int nuNphiBasis = 0;
	int phiNnuBasis = 0;
	int phiNphiBasis = 0;

	fourierDef0.setSize( refAdjust_muinNnuBasis, refAdjust_muinNphiBasis, nuNnuBasis,
		 nuNphiBasis, phiNnuBasis, phiNphiBasis );

	delta_a = 0;
}


void LVReferenceShape::setupRegularRefSurfaces( vector<double> & cin0In, vector <double> & cout0In,
		int NnuSplines, int NphiSplines, double nuUp0MinIn, double muSplineConst, double aIn )
{

	nuUp0Min = nuUp0MinIn;

	// Copy the vector of parameters
	cin0.resize( cin0In.size() );
	cout0.resize( cout0In.size() );
	for(int i = 0; i < cin0.size(); i++)
	{
		cin0[i] = cin0In[i];
	}
	for(int i = 0; i < cout0.size(); i++)
	{
		cout0[i] = cout0In[i];
	}

	// Now compute the surfaces
	endo0.createRegularSplines( cin0.data(), NnuSplines, NphiSplines, nuUp0Min, muSplineConst );
	epi0.createRegularSplines( cout0.data(), NnuSplines, NphiSplines, nuUp0Min, muSplineConst );

	// Normally this is the same
	a_original = aIn;

}
