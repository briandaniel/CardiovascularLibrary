/*
 * LVSplineDeformation.cpp
 *
 *  Created on: Jun 5, 2018
 *      Author: brian
 */

#include "LVBicubicDeformation.hpp"





// Main function
// first are the coordinate inputs
// second are the precomputed parameters
// third, the outputs required by the regular program
void BicubicDeformation::computeDeformation( double muin0, double nu0, double phi0,
	double nuUp0, double fc0, double acube, double chmu0, double shmu0, double chmu0sq, double cnu0sq, double cphi0, double sphi0,
	double snu0, double snu0sq, double cnu0, double dfc0dnu0, double dfc0dphi0, double dmuin0dnu0, double dmuin0dphi0,
	double &mu, double &nu, double &phi, double &shmu, double &chmu, double &shmusq, double &chmusq, double &snu, double &cnu, double &snusq, double &cnusq,
	double &fc, double &dfcdnu0, double &dfcdphi0, double &K, double &sinRatio, double & dnudnu0, double &dnudphi0, double &dphidnu0, double &dphidphi0,
	double &d2nudnu02, double &d2nudphi0dnu0, double & d2nudphi02, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0,
	double &dmudmu0, double &dmudnu0, double & dmudphi0 )
{

	// local variables
	double R, Knu0, Kphi0, A, Bnu0, Bphi0, Rnu0, Rphi0;

	//--------------------------------------------------------------------//
	// 1. compute nu/phi and derivatives
	//--------------------------------------------------------------------//

	// Compute nu and nu derivatives
	nu = nuFuncBicub( nu0, phi0, nuSpline );
	nuDerivativesBicub( nu0, phi0, nuSpline, dnudnu0, dnudphi0 );
	nuSecondDerivativesBicub( nu0, phi0, nuSpline, d2nudnu02, d2nudphi0dnu0, d2nudphi02 );

	// Compute phi and phi derivatives
	phi = phiFuncBicub( nu0, phi0, phiSpline );
	phiDerivativesBicub( nu0, phi0, phiSpline, dphidnu0, dphidphi0 );
	phiSecondDerivativesBicub( nu0, phi0, phiSpline, d2phidnu02, d2phidphi02, d2phidphi0dnu0 );

	//--------------------------------------------------------------------//
	// 2. compute mu
	//--------------------------------------------------------------------//

	// Need a function for this because of the singularity at nu = pi
	sinRatio = computeSinRatioBicub( nu0, nu );

	// compute K, R
	K = dnudnu0*dphidphi0 - dnudphi0*dphidnu0;
	R = sinRatio*K;

	// compute functions of nu/phi
	snu = sin(nu);
	cnu = cos(nu);
	snusq = pow(snu,2);
	cnusq = pow(cnu,2);

	// The contraction function
	fc = fcfuncBicub( muin0, nu0, phi0, fmuinSpline, fc0, acube, R, cnusq);

	// Compute mu
	mu = muFuncBicub( K, sinRatio, nu, acube, chmu0, chmu0sq, cnu0sq, fc  );

	// compute functions of mu
	shmu = sinh(mu);
	chmu = cosh(mu);
	shmusq = pow(shmu,2);
	chmusq = pow(chmu,2);


	//--------------------------------------------------------------------//
	// 3. compute mu derivatives
	//--------------------------------------------------------------------//

	// aux variables
	Knu0 = d2nudnu02*dphidphi0 + dnudnu0*d2phidphi0dnu0 - d2nudphi0dnu0*dphidnu0 - dnudphi0*d2phidnu02;
	Kphi0 = d2nudphi0dnu0*dphidphi0 + dnudnu0*d2phidphi02 - d2nudphi02*dphidnu0 - dnudphi0*d2phidphi0dnu0;

	Rnu0 = ( snu0*cnu*dnudnu0 - snu*cnu0 )/( snu0sq )*K + sinRatio*Knu0;
	Rphi0 = ( cnu/snu0 )*dnudphi0*K + sinRatio*Kphi0;

	// derivatives of fc
	fcDerivativesBicub( muin0, nu0, phi0, fmuinSpline, dfc0dnu0, dfc0dphi0, acube,
			dmuin0dnu0, dmuin0dphi0, R, Rnu0, Rphi0, cnu, cnusq, snu, dnudnu0,dnudphi0,
			dfcdnu0, dfcdphi0 );

	// more aux variables
	A = THIRD*pow(chmu,3) - chmu*cnusq - THIRD + cnusq;
	Bnu0 = 2*cnu0*snu0*(chmu0 - 1.0) - 1.0/acube*dfcdnu0;
	Bphi0 = -1.0/acube*dfcdphi0;

	// mu derivatives
	dmudmu0 = (1.0/sinRatio)*( shmu0*( chmu0sq - cnu0sq ) )/( K*shmu*( chmusq - cnusq ) );

	dmudnu0 = ( Bnu0 + 2*R*cnu*snu*dnudnu0*(1-chmu) - Rnu0*A )/( R*shmu*(chmusq - cnusq) );

	dmudphi0 = ( Bphi0 + 2*R*cnu*snu*dnudphi0*(1-chmu) - Rphi0*A )/( R*shmu*(chmusq - cnusq) );


	// correction near nu0 = pi due to singularity of dmudnu0
	// assumes construction with dmudnu0 = 0 at nu = pi
	if( fabs(nu0 - PI) < 1e-6)
	{
		dmudnu0 = 0;
	}

}





// Pulls the values of q
void pullq_parameters_bicub( double *q, int Nc, int Nd, int Ne,  double *sc, double *sd, double *se)
{

	//
	sc[0] = q[0];
	sc[1] = 0.0;
	sc[2] = 0.0;
	for (int k = 1; k < Nc; k++ )
		sc[k+2] = q[k];



	// set the first three values to zero, and then fill in the remaining values
	// because delta phi and derivatives are zero at the apex
	sd[0] = 0.0;
	sd[1] = 0.0;
	sd[2] = 0.0;
	for (int k = 0; k < Nd; k++ )
		sd[k+3] = q[k+Nc];


	// set the first three values to zero, and then fill in the remaining values
	// because delta nu and derivatives are zero at the apex
	se[0] = 0.0;
	se[1] = 0.0;
	se[2] = 0.0;
	for (int k = 0; k < Ne; k++ )
		se[k+3] = q[k+Nc+Nd];



}

void BicubicDeformation::updateBicubicDeformationSplines( double * q, int Nc, int Nd, int Ne )
{

	// pull the variables
	double *sc = new double[ Nc + 2];
	double *sd = new double[ Nd + 3];
	double *se = new double[ Ne + 3];
	pullq_parameters_bicub( q, Nc, Nd, Ne, sc, sd, se );

	// update the spline parameters
	fmuinSpline.computeSplineParameters( sc );
	phiSpline.computeSplineParameters( sd );
	nuSpline.computeSplineParameters( se );

	// clean up
	delete [] sc;
	delete [] sd;
	delete [] se;

}



void BicubicDeformation::setupSplines( int muinNnuGrid, int muinNphiGrid, int nuNnuGrid, int nuNphiGrid,
		int phiNnuGrid, int phiNphiGrid, double nuUp0Min, double muin0_const )
{

	// create splines with the correct sizes
	fmuinSpline.createConstantProlateSplines(0.0, muinNnuGrid, muinNphiGrid, nuUp0Min, muin0_const );

	nuSpline.createConstantProlateSplines(0.0, nuNnuGrid, nuNphiGrid, nuUp0Min, muin0_const );

	phiSpline.createConstantProlateSplines(0.0, phiNnuGrid, phiNphiGrid, nuUp0Min, muin0_const );


}


BicubicDeformation::BicubicDeformation()
{

}




BicubicDeformation::BicubicDeformation( int muinNnuGrid, int muinNphiGrid, int nuNnuGrid, int nuNphiGrid,
		int phiNnuGrid, int phiNphiGrid, double nuUp0Min, double muin0_const )
{

	// create splines with the correct sizes
	fmuinSpline.createConstantProlateSplines(0.0, muinNnuGrid, muinNphiGrid, nuUp0Min, muin0_const );

	nuSpline.createConstantProlateSplines(0.0, nuNnuGrid, nuNphiGrid, nuUp0Min, muin0_const );

	phiSpline.createConstantProlateSplines(0.0, phiNnuGrid, phiNphiGrid, nuUp0Min, muin0_const );


}

BicubicDeformation::~BicubicDeformation()
{
}


// nu functions
double nuFuncBicub( double nu0, double phi0, ProlateSplines & nuSpline)
{

	double nu = nu0 + nuSpline.evaluateSplines(nu0,phi0);
	return nu;
}

void nuDerivativesBicub( double nu0, double phi0, ProlateSplines & nuSpline, double &dnudnu0, double &dnudphi0 )
{

	nuSpline.evaluateSplineDerivatives(nu0,phi0,dnudnu0,dnudphi0);
	dnudnu0 = 1 + dnudnu0;

}

void nuSecondDerivativesBicub( double nu0, double phi0, ProlateSplines & nuSpline, double &d2nudnu02, double &d2nudphi0dnu0, double &d2nudphi02 )
{

	nuSpline.evaluateSplineSecondDerivatives( nu0, phi0, d2nudnu02, d2nudphi0dnu0, d2nudphi02 );

}

double computeSinRatioBicub( double nu0, double nu ){

	double sinRatio = sin(nu)/sin(nu0);

	// The sin ratio is singular at nu0 = pi, but by the definitions of
	// the transformation the ratio is 1
	if( fabs(nu0 - PI) < 1e-7)
	{
		sinRatio = 1;
	}

	return sinRatio;
}


// phi functions
double phiFuncBicub(  double nu0, double phi0, ProlateSplines & phiSpline)
{
	double phi = phi0 + phiSpline.evaluateSplines(nu0,phi0);
	return phi;
}

void phiDerivativesBicub( double nu0, double phi0, ProlateSplines & phiSpline, double &dphidnu0, double &dphidphi0 )
{
	phiSpline.evaluateSplineDerivatives(nu0,phi0,dphidnu0,dphidphi0);
	dphidphi0 = 1 + dphidphi0;

}

void phiSecondDerivativesBicub( double nu0, double phi0, ProlateSplines & phiSpline, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0 )
{
	phiSpline.evaluateSplineSecondDerivatives( nu0, phi0, d2phidnu02, d2phidphi0dnu0, d2phidphi02 );

}




double fcfuncBicub( double muin0, double nu0, double phi0, ProlateSplines & fmuinSpline, double fc0, double acube, double R, double cnusq)
{
	double fmuin, muin, chmuin, fc;

	fmuin = fmuinSpline.evaluateSplines(nu0,phi0);
	muin = muin0 + fmuin;
	chmuin = cosh(muin);

	fc = fc0 - acube*R*( chmuin *(THIRD * pow(chmuin,2) - cnusq) - (THIRD - cnusq) );

	// cant contract past fc0, which is contraction down to a single line
	if ( fc > fc0 ){ fc = fc0; }

	return fc;

}



double muFuncBicub( double K, double sinRatio, double nu, double acube, double chmu0, double chmu0sq, double cnu0sq, double fc  )
{

	double b, R, numerator, S;
	double mu;

	std::complex <double> root3, coshMuComplex;


	R = K*sinRatio;
	b = - pow( cos(nu), 2);
	numerator = acube*( chmu0 * ( THIRD * chmu0sq - cnu0sq ) - ( THIRD - cnu0sq) ) - fc;
	S = - numerator/(acube*R) - ( THIRD - pow(cos(nu),2) );

	// This has to be computed with complex numbers
	root3 = pow( sqrt( (std::complex<double>)  ( 4*pow(b,3) + 9*pow(S,2) ) ) - 3*S, 1.0/3.0 );
	coshMuComplex = root3/RTHIRD - RTHIRD*b/root3;

	mu = acosh( coshMuComplex.real() );

	// in case of over contraction simply set mu = 0
	if(mu != mu)
		mu = 0.0;

	return mu;
}




void fcDerivativesBicub(double muin0, double nu0, double phi0, ProlateSplines & fmuinSpline, double dfc0dnu0, double dfc0dphi0, double acube,
		double dmuin0dnu0, double dmuin0dphi0, double R, double Rnu0, double Rphi0, double cnu, double cnusq, double snu, double dnudnu0, double dnudphi0,
		double & dfcdnu0, double & dfcdphi0 )
{

	double fmuin, dfmuindnu0, dfmuindphi0, muin, dmuindnu0, dmuindphi0, W, dWdnu0, dWdphi0;
	double chmuin, chmuinsq, shmuin;

	fmuin = fmuinSpline.evaluateSplines(nu0,phi0);
	muin = muin0 + fmuin;
	chmuin = cosh(muin);
	chmuinsq = pow(chmuin,2);
	shmuin = sinh(muin);


	// muin derivatives
	fmuinSpline.evaluateSplineDerivatives(nu0, phi0, dfmuindnu0, dfmuindphi0);
	dmuindnu0 = dmuin0dnu0 + dfmuindnu0;
	dmuindphi0 = dmuin0dphi0 + dfmuindphi0;


	// aux variable W and derivatives
	W = chmuin*( THIRD*chmuinsq - cnusq ) - (THIRD - cnusq);
	dWdnu0 = shmuin*dmuindnu0*( chmuinsq - cnusq ) + 2*cnu*snu*dnudnu0*( chmuin - 1 );
	dWdphi0 = shmuin*dmuindphi0*( chmuinsq - cnusq ) + 2*cnu*snu*dnudphi0*( chmuin - 1 );


	// final fc derivatives
	dfcdnu0 = dfc0dnu0 - acube*( Rnu0*W + R*dWdnu0 );
	dfcdphi0 = dfc0dphi0 - acube*( Rphi0*W + R*dWdphi0 );

}











