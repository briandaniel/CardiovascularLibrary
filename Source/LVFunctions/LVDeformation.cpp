/*
 * deformation.cpp
 *
 *  Created on: Sep 12, 2017
 *      Author: brian
 */

#include "../LVFunctions/LVDeformation.hpp"

// transition location of the nu cutoffs for phi modes; easy to define it here
#define NU0T 2.3562
#define NUR 2.7489
#define LVAL -1.5708




// Main function
// first are the coordinate inputs
// second are the precomputed parameters
// third, the outputs required by the regular program
void computeSevenParamDeformation( double muin0, double nu0, double phi0, double * q,
		double nuUp0, double fc0, double acube, double chmu0, double shmu0, double chmu0sq, double cnu0sq, double cphi0, double sphi0,
		double snu0, double snu0sq, double cnu0, double dfc0dnu0, double dfc0dphi0, double dmuin0dnu0, double dmuin0dphi0, double chmue,
		double &mu, double &nu, double &phi, double &shmu, double &chmu, double &shmusq, double &chmusq, double &snu, double &cnu, double &snusq, double &cnusq,
		double &fc, double &dfcdnu0, double &dfcdphi0, double &K, double &sinRatio, double & dnudnu0, double &dnudphi0, double &dphidnu0, double &dphidphi0,
		double &d2nudnu02, double &d2nudphi0dnu0, double & d2nudphi02, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0,
		double &dmudmu0, double &dmudnu0, double & dmudphi0)
{

	// local variables
	int Nc, Nd, Ne, Nq;
	double R, Knu0, Kphi0, A, Bnu0, Bphi0, Rnu0, Rphi0;

	// aux variables
	currentDeformationParameterSize(Nc, Nd, Ne, Nq );
	double *c = new double[ Nc ];
	double *d = new double[ Nd ];
	double *e = new double[ Ne ];
	pullq_parameters( q, Nc, Nd, Ne, c, d, e );


	//--------------------------------------------------------------------//
	// 1. compute nu/phi and derivatives
	//--------------------------------------------------------------------//

	// Compute nu and nu derivatives
	nu = nuFunc( nu0, phi0, nuUp0, e );
	nuDerivatives( nu0, phi0, nuUp0, e, dnudnu0, dnudphi0 );
	nuSecondDerivatives( nu0, phi0, nuUp0, e, d2nudnu02, d2nudphi0dnu0, d2nudphi02 );

	// Compute phi and phi derivatives
	phi = phiFunc(  nu0, phi0, nuUp0,  d );
	phiDerivatives( nu0, phi0, nuUp0, d, dphidnu0, dphidphi0 );
	phiSecondDerivatives( nu0, phi0, nuUp0, d, d2phidnu02, d2phidphi02, d2phidphi0dnu0 );


	//--------------------------------------------------------------------//
	// 2. compute mu
	//--------------------------------------------------------------------//

	// Need a function for this because of the singularity at nu = pi
	sinRatio = computeSinRatio( nu0, nu );

	// compute K, R
	K = dnudnu0*dphidphi0 - dnudphi0*dphidnu0;
	R = sinRatio*K;

	// compute functions of nu/phi
	snu = sin(nu);
	cnu = cos(nu);
	snusq = pow(snu,2);
	cnusq = pow(cnu,2);

	// The contraction function
	fc = ffunc(  c, acube, cphi0, sphi0, nu0, fc0, chmue, K, sinRatio, nu );

	// Compute mu
	mu = muFunc( K, sinRatio, nu, acube, chmu0, chmu0sq, cnu0sq, fc  );

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
	fDerivatives(  acube,  sphi0,  cphi0,  nu0, fc0,  dfc0dnu0,  dfc0dphi0,  R,  Rphi0,
			 Rnu0,  chmue,  snu,  cnu,  dnudphi0,  dnudnu0,  c,  dfcdnu0,  dfcdphi0 );

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



	// clean up
	delete [] c;
	delete [] d;
	delete [] e;
}



// Gives the size of the vectors in the current deformation computation
void currentDeformationParameterSize( int & Nc, int & Nd, int & Ne, int & Nq )
{
	Nc = 4;
	Nd = 1;
	Ne = 2;

	Nq = Nc + Nd + Ne;
}




// Pulls the values of q
void pullq_parameters( double *q, int Nc, int Nd, int Ne,  double *c, double *d, double *e)
{

	for (int k = 0; k < Nc; k++ )
		c[k] = q[k];

	for (int k = 0; k < Nd; k++ )
		d[k] = q[k+Nc];

	for (int k = 0; k < Ne; k++ )
		e[k] = q[k+Nc+Nd];

}



double nuFunc( double nu0, double phi0, double nuUp0, double * e )
{

	double v1, v2, mv1, mv2, nu1;
	double nu;

	v1 = pow(nu0 - PI,2);
	v2 = pow(nu0 - PI,3);

	mv1 = ( 12*e[0] + 2*PI*e[1] )/(pow(PI,2));
	mv2 = ( 16*e[0] + 4*PI*e[1] )/pow(PI,3);

	nu1 = mv1*v1 + mv2*v2;
	nu = nu0 + nu1;

	return nu;
}


void nuDerivatives( double nu0, double phi0, double nuUp0, double * e, double &dnudnu0, double &dnudphi0 )
{

	double v1dnu0, v2dnu0, mv1, mv2;
	double nu1_dnu0;


	v1dnu0 = 2*(nu0 - PI);
	v2dnu0 = 3*pow(nu0 - PI,2);

	mv1 = ( 12*e[0] + 2*PI*e[1] )/(pow(PI,2));
	mv2 = ( 16*e[0] + 4*PI*e[1] )/pow(PI,3);

	nu1_dnu0 = mv1*v1dnu0 + mv2*v2dnu0;

	dnudnu0 = 1 + nu1_dnu0;
	dnudphi0 = 0;

}



void nuSecondDerivatives( double nu0, double phi0, double nuUp0, double * e, double &d2nudnu02, double &d2nudphi0dnu0, double &d2nudphi02 )
{

	double v1dnu0, v2dnu0, v1dnu02, v2dnu02, mv1, mv2, d2nu1_dnu02;


	v1dnu0 = 2*(nu0 - PI);
	v2dnu0 = 3*pow(nu0 - PI,2);
	v1dnu02 = 2;
	v2dnu02 = 6*(nu0 - PI);

	mv1 = ( 12*e[0] + 2*PI*e[1] )/(pow(PI,2));
	mv2 = ( 16*e[0] + 4*PI*e[1] )/pow(PI,3);

	d2nu1_dnu02 = mv1*v1dnu02 + mv2*v2dnu02;

	d2nudnu02 = d2nu1_dnu02;
	d2nudphi0dnu0 = 0;
	d2nudphi02 = 0;

}



double computeSinRatio( double nu0, double nu ){

	double sinRatio = sin(nu)/sin(nu0);

	// The sin ratio is singular at nu0 = pi, but by the definitions of
	// the transformation the ratio is 1
	if( fabs(nu0 - PI) < 1e-7)
	{
		sinRatio = 1;
	}

	return sinRatio;
}



double phiFunc(  double nu0, double phi0, double nuUp0, double * d )
{
	double nur, deltaPhi;
	double phi;

	nur = NUR;

	if( nu0 <= nur )
	{
		deltaPhi = d[0]*(nu0 - nur);
	}else
	{
		deltaPhi = d[0]/( 2*(nur-PI) )*pow(nu0-PI,2) - d[0]*(nur-PI)/2;
	}

	phi = phi0+deltaPhi;
	return phi;
}



void phiDerivatives( double nu0, double phi0, double nuUp0, double * d, double &dphidnu0, double &dphidphi0 )
{
	double nur, deltaPhi_dnu0, deltaPhi_dphi0;

	nur = NUR;

	if( nu0 < nur ){
		deltaPhi_dnu0 = d[0];
	}else{
		deltaPhi_dnu0 = 2*d[0]/( 2*(nur-PI) )*(nu0-PI);
	}

	deltaPhi_dphi0 = 0;

	dphidnu0 = deltaPhi_dnu0;
	dphidphi0 = 1 + deltaPhi_dphi0;

}


void phiSecondDerivatives( double nu0, double phi0, double nuUp0, double * d, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0 )
{
	double nur, deltaPhi_dnu02;

	nur = NUR;

	if( nu0 < nur ){
		deltaPhi_dnu02 = 0;
	}else{
		deltaPhi_dnu02 = 2*d[0]/( 2*(nur-PI) );
	}

	d2phidnu02 = deltaPhi_dnu02;
	d2phidphi02 = 0;
	d2phidphi0dnu0 = 0;

}



double ffunc( double * c, double acube, double cphi0, double sphi0, double nu0,
						  double fc0, double chmue, double K, double sinRatio, double nu)
{
	double fseries_phi0, Gnu0, fseries, fc1, fc2, fe, cnusq;
	double fc;


	// Modes for the regular function
	Gnu0 = 0.5*( 1 + cos( PI*(nu0 - NU0T)/(NU0T-PI) ) )*(nu0 >= NU0T) + ( nu0 < NU0T  );
	fc1 = fc0*(c[0] + c[1]*Gnu0*cphi0 + c[2]*Gnu0*sphi0);

	// this adjusts toward the elliptical shape mue
	cnusq = pow( cos(nu), 2);
	fc2 = c[3]*( fc0 - acube*K*sinRatio*( chmue*( THIRD*pow(chmue,2) - cnusq  ) - ( THIRD - cnusq ) ) );

	// sum
	fc = fc1 + fc2;

	// cant contract past f0, which is contraction down to a single line
	if ( fc > fc0 ){ fc = fc0; }

	return fc;

}


void fDerivatives( double acube, double sphi0, double cphi0, double nu0,double fc0, double dfc0dnu0, double dfc0dphi0, double R, double Rphi0,
		double Rnu0, double chmue, double snu, double cnu, double dnudphi0, double dnudnu0, double * c, double &dfcdnu0, double &dfcdphi0 )
{

	double cnusq;
	double Gnu0, dGnu0dnu0, fc1;
	double dfcCoef, dfcCoef_dnu0, dfcCoef_dphi0;
	double dfseries_phi0_dphi0, dfseries_nu0_dnu0;
	double fseries, dfseries_dnu0, dfseries_dphi0;
	double feRless, dfeRless_dnu0, dfeRless_dphi0;
	double dfe_dnu0, dfe_dphi0;
	double dfc1_dnu0, dfc1_dphi0;
	double dfc2_dnu0, dfc2_dphi0;

	// fseries derivatives -- must match the values in ffunc above

	// fc1 derivative
	Gnu0 = 0.5*( 1 + cos( PI*(nu0 - NU0T)/(NU0T-PI) ) )*(nu0 >= NU0T) + ( nu0 < NU0T  );
	dGnu0dnu0 = PI/(2*(PI-NU0T)) * ( sin(PI*(nu0 - NU0T)/(NU0T-PI)  ) )*(nu0 >= NU0T);

	dfcCoef = (c[0] + c[1]*Gnu0*cphi0 + c[2]*Gnu0*sphi0);
	dfcCoef_dnu0 = c[1]*dGnu0dnu0*cphi0 + c[2]*dGnu0dnu0*sphi0;
	dfcCoef_dphi0 = -c[1]*Gnu0*sphi0 + c[2]*Gnu0*cphi0;

	dfc1_dnu0 = dfc0dnu0*dfcCoef + fc0*dfcCoef_dnu0;
	dfc1_dphi0 = dfc0dphi0*dfcCoef + fc0*dfcCoef_dphi0;


	// fc2 derivative

	cnusq = pow( cnu, 2);
	feRless = acube*( chmue*( THIRD*pow(chmue,2) - cnusq  ) - ( THIRD - cnusq ) ) ;
	dfeRless_dnu0 = 2*acube*snu*cnu*dnudnu0*(chmue - 1);
	dfeRless_dphi0 = 2*acube*snu*cnu*dnudphi0*(chmue - 1);
	dfe_dnu0 = Rnu0*feRless + R*dfeRless_dnu0;
	dfe_dphi0 = Rphi0*feRless + R*dfeRless_dphi0;

	dfc2_dnu0 = c[3]*( dfc0dnu0 - dfe_dnu0);
	dfc2_dphi0 = c[3]*( dfc0dphi0 - dfe_dphi0);



	// sum
	dfcdnu0 = dfc1_dnu0 + dfc2_dnu0;
	dfcdphi0 = dfc1_dphi0 + dfc2_dphi0;


}



double muFunc( double K, double sinRatio, double nu, double acube, double chmu0, double chmu0sq, double cnu0sq, double fc  )
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






















