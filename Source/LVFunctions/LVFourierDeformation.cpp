/*
 * LVFourierDeformation.cpp
 *
 *  Created on: Jun 17, 2018
 *      Author: brian
 */

#include "LVFourierDeformation.hpp"



/*
 * THIS IS THE ORIGINAL VERSION
 *
*/
// Main function
// first are the coordinate inputs
// second are the precomputed parameters
// third, the outputs required by the regular program
void FourierDeformation::computeDeformation( double muin0, double nu0, double phi0, double * q,
	double nuUp0, double fc0, double acube, double chmu0, double shmu0, double chmu0sq, double cnu0sq, double cphi0, double sphi0,
	double snu0, double snu0sq, double cnu0, double dfc0dnu0, double dfc0dphi0, double dmuin0dnu0, double dmuin0dphi0,
	double &mu, double &nu, double &phi, double &shmu, double &chmu, double &shmusq, double &chmusq, double &snu, double &cnu, double &snusq, double &cnusq,
	double &fc, double &dfcdnu0, double &dfcdphi0, double &K, double &sinRatio, double & dnudnu0, double &dnudphi0, double &dphidnu0, double &dphidphi0,
	double &d2nudnu02, double &d2nudphi0dnu0, double & d2nudphi02, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0,
	double &dmudmu0, double &dmudnu0, double & dmudphi0 )
{



	double muin, dmuindnu0, dmuindphi0;

	// These are almost the same as c, d, and e
	// but have one extra zero for d, e
	double *coef_cf = new double[ Nc ];
	double *coef_df = new double[ Nd + 1];
	double *coef_ef = new double[ Ne + 1];
	pullq_parameters( q, coef_cf, coef_df, coef_ef );


	// local variables
	double R, Knu0, Kphi0, A, Bnu0, Bphi0, Rnu0, Rphi0;

	//--------------------------------------------------------------------//
	// 1. compute nu/phi and derivatives
	//--------------------------------------------------------------------//

	// Compute nu and nu derivatives
	nuFunc( nu0, phi0, coef_ef, nuUp0, nu, dnudnu0,  dnudphi0,  d2nudnu02,  d2nudphi0dnu0,  d2nudphi02 );

	// Compute phi and phi derivatives
	phiFunc( nu0, phi0, coef_df, nuUp0, phi, dphidnu0, dphidphi0, d2phidnu02, d2phidphi0dnu0, d2phidphi02 );

	//--------------------------------------------------------------------//
	// 2. compute mu
	//--------------------------------------------------------------------//

	// Need a function for this because of the singularity at nu = pi
	sinRatio = computeSinRatioFourier( nu0, nu );

	// compute K, R
	K = dnudnu0*dphidphi0 - dnudphi0*dphidnu0;
	R = sinRatio*K;

	// compute functions of nu/phi
	snu = sin(nu);
	cnu = cos(nu);
	snusq = pow(snu,2);
	cnusq = pow(cnu,2);

	// The contraction function
	muinFunc( muin0, nu0, phi0, coef_cf, dmuin0dnu0, dmuin0dphi0,  nuUp0, muin, dmuindnu0, dmuindphi0 );
	fc = fcfuncFourier( muin, fc0, acube, R, cnusq);

	// Compute mu
	mu = muFuncFourier( K, sinRatio, nu, acube, chmu0, chmu0sq, cnu0sq, fc  );


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
	fcDerivativesFourier( muin, dmuindnu0, dmuindphi0, dfc0dnu0, dfc0dphi0, acube,
			R,  Rnu0, Rphi0,cnu, cnusq, snu, dnudnu0, dnudphi0, dfcdnu0, dfcdphi0 );

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


	delete [] coef_cf;
	delete [] coef_df;
	delete [] coef_ef;

	//		cout << "nu0 = " << nu0 << ", phi0 =" << phi0 << ", mu0 = " << asinh(shmu0) << ", mu = " << mu << ", fc = " << fc << ", fc0 = " << fc0 << endl;
}




void computeFourierValueAndDerivatives( double nu0, double phi0, double * coef, int Nnu0, int Nphi0, double nuUp0,
		double & f, double & dfdnu0, double & dfdphi0, double &d2fdnu02, double &d2fdphi0dnu0, double &d2fdphi02 )
{


	double fnu0, fcphi0, fsphi0;
	double dfnu0_dnu0, dfcphi0_dphi0, dfsphi0_dphi0;
	double d2fnu0_dnu02, d2fcphi0_dphi02, d2fsphi0_dphi02;

	int k = 0;
	// constant value first
	f = coef[k];
	k++;

	dfdnu0 = 0;
	dfdphi0 = 0;

	d2fdnu02 = 0;
	d2fdphi0dnu0 = 0;
	d2fdphi02 = 0;

	// 1. Constant phi0 modes with variations in nu0
	for (int i = 0; i < Nnu0; i++)
	{
		int jcos = 0;
		int jsin = 0;
		computeFourierFunctionValuesAndDerivativesSinglePoint( nu0, phi0, i, jcos, jsin, nuUp0,
				fnu0, dfnu0_dnu0, d2fnu0_dnu02, fcphi0, dfcphi0_dphi0, d2fcphi0_dphi02,	fsphi0, dfsphi0_dphi0, d2fsphi0_dphi02 );

		f = f + coef[k]*fnu0*fcphi0;

		dfdnu0 = dfdnu0 + coef[k]*dfnu0_dnu0*fcphi0;
		dfdphi0 = dfdphi0 + coef[k]*fnu0*dfcphi0_dphi0;

		d2fdnu02 = d2fdnu02 + coef[k]*d2fnu0_dnu02*fcphi0;
		d2fdphi0dnu0 = d2fdphi0dnu0 + coef[k]*dfnu0_dnu0*dfcphi0_dphi0;
		d2fdphi02 = d2fdphi02 + coef[k]*fnu0*d2fcphi0_dphi02;

		k++;
	}


	// 2. Cycle through cos/sin modes
	for (int j = 0; j < Nphi0; j++)
	{
		for (int i = 0; i < Nnu0; i++)
		{

			int jcos = j+1;
			int jsin = j+1;
			computeFourierFunctionValuesAndDerivativesSinglePoint( nu0, phi0, i, jcos, jsin, nuUp0,
					fnu0, dfnu0_dnu0, d2fnu0_dnu02, fcphi0, dfcphi0_dphi0, d2fcphi0_dphi02,	fsphi0, dfsphi0_dphi0, d2fsphi0_dphi02 );

		    // cosine mode
			f = f + coef[k]*fnu0*fcphi0;

			dfdnu0 = dfdnu0 + coef[k]*dfnu0_dnu0*fcphi0;
			dfdphi0 = dfdphi0 + coef[k]*fnu0*dfcphi0_dphi0;

			d2fdnu02 = d2fdnu02 + coef[k]*d2fnu0_dnu02*fcphi0;
			d2fdphi0dnu0 = d2fdphi0dnu0 + coef[k]*dfnu0_dnu0*dfcphi0_dphi0;
			d2fdphi02 = d2fdphi02 + coef[k]*fnu0*d2fcphi0_dphi02;
			k++;

			// sine mode
			f = f + coef[k]*fnu0*fsphi0;

			dfdnu0 = dfdnu0 + coef[k]*dfnu0_dnu0*fsphi0;
			dfdphi0 = dfdphi0 + coef[k]*fnu0*dfsphi0_dphi0;

			d2fdnu02 = d2fdnu02 + coef[k]*d2fnu0_dnu02*fsphi0;
			d2fdphi0dnu0 = d2fdphi0dnu0 + coef[k]*dfnu0_dnu0*dfsphi0_dphi0;
			d2fdphi02 = d2fdphi02 + coef[k]*fnu0*d2fsphi0_dphi02;
			k++;


		}
	}


}


void computeFourierFunctionValuesAndDerivativesSinglePoint( double nu0, double phi0, int i, int jcos, int jsin, double nuUp0,
		double & fnu0, double &  dfnu0_dnu0, double & d2fnu0_dnu02,
		double & fcphi0, double & dfcphi0_dphi0, double & d2fcphi0_dphi02,
		double & fsphi0, double & dfsphi0_dphi0, double & d2fsphi0_dphi02 )
{

	fnu0 = 1 - cos( (i+1)*(nu0 - PI));
	dfnu0_dnu0 = (i+1)*sin( (i+1)*(nu0 - PI));
	d2fnu0_dnu02 = (i+1)*(i+1)*cos( (i+1)*(nu0 - PI) );

	fcphi0 = cos(jcos*phi0);
	dfcphi0_dphi0 = -sin(jcos*phi0)*jcos;
	d2fcphi0_dphi02 = -cos(jcos*phi0)*jcos*jcos;

	fsphi0 = sin(jsin*phi0);
	dfsphi0_dphi0 = cos(jsin*phi0)*jsin;
	d2fsphi0_dphi02 = -sin(jsin*phi0)*jsin*jsin;

}


// Pulls the values of q
void FourierDeformation::pullq_parameters( double *q, double *coef_cf, double *coef_df, double *coef_ef )
{

	//
	for (int k = 0; k < Nc; k++ )
	{
		coef_cf[k] = q[k];
	}

	// set the first value to zero, and then fill in the remaining values
	// because delta phi is zero at the apex
	coef_df[0] = 0.0;
	for (int k = 0; k < Nd; k++ )
		coef_df[k+1] = q[k+Nc];

	// set the first value to zero, and then fill in the remaining values
	// because delta nu is zero at the apex
	coef_ef[0] = 0.0;
	for (int k = 0; k < Ne; k++ )
		coef_ef[k+1] = q[k+ Nc + Nd];

}



void FourierDeformation::setSize( int muinNnuBasisIn, int muinNphiBasisIn, int nuNnuBasisIn,
		int nuNphiBasisIn, int phiNnuBasisIn, int phiNphiBasisIn )
{

	muinNnuBasis = muinNnuBasisIn;
	muinNphiBasis = muinNphiBasisIn;
	Nc = muinNnuBasis*( muinNphiBasis*2+1 ) + 1;

	nuNnuBasis = nuNnuBasisIn;
	nuNphiBasis = nuNphiBasisIn;
	Ne = nuNnuBasis*( nuNphiBasis*2+1 );

	phiNnuBasis = phiNnuBasisIn;
	phiNphiBasis = phiNphiBasisIn;
	Nd = phiNnuBasis*( phiNphiBasis*2+1 );



}


void FourierDeformation::computeNParams( int & NqOut, int & NcOut, int & NdOut, int & NeOut )
{

	NcOut = Nc;
	NdOut = Nd;
	NeOut = Ne;

	NqOut = Nc+Nd+Ne;

}

FourierDeformation::FourierDeformation()
{
	setSize( 1, 1, 1, 1, 1, 1 );
}

void FourierDeformation::printNModesToConsole()
{
	cout << "Nc, Nd, Ne = " << Nc << ", " << Nd << ", " << Ne << endl;

	cout << " muinNnuBasis = " << muinNnuBasis << endl;
	cout << " muinNphiBasis = " << muinNphiBasis << endl;

	cout << " nuNnuBasis = " << nuNnuBasis << endl;
	cout << " nuNphiBasis = " << nuNphiBasis << endl;

	cout << " phiNnuBasis = " << phiNnuBasis << endl;
	cout << " phiNphiBasis = " << phiNphiBasis << endl;

}





FourierDeformation::FourierDeformation( int muinNnuBasisIn, int muinNphiBasisIn, int nuNnuBasisIn,
		int nuNphiBasisIn, int phiNnuBasisIn, int phiNphiBasisIn )
{

	setSize( muinNnuBasisIn, muinNphiBasisIn, nuNnuBasisIn,
		 nuNphiBasisIn, phiNnuBasisIn, phiNphiBasisIn );

}


void FourierDeformation::nuFunc( double nu0, double phi0, double * ef, double nuUp0,
			double & nu, double & dnudnu0, double & dnudphi0, double & d2nudnu02, double & d2nudphi0dnu0, double & d2nudphi02 )
{

	double f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02;

	computeFourierValueAndDerivatives( nu0, phi0, ef, nuNnuBasis, nuNphiBasis, nuUp0,
			f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02 );

	nu = nu0 + f;
	dnudnu0 = 1 + dfdnu0;
	dnudphi0 = dfdphi0;
	d2nudnu02 = d2fdnu02;
	d2nudphi0dnu0 = d2fdphi0dnu0;
	d2nudphi02 = d2fdphi02;

}



double computeSinRatioFourier( double nu0, double nu ){

	double sinRatio = sin(nu)/sin(nu0);

	// The sin ratio is singular at nu0 = pi, but by the definitions of
	// the transformation the ratio is 1
	if( fabs(nu0 - PI) < 1e-7)
	{
		sinRatio = 1;
	}

	return sinRatio;
}



void FourierDeformation::phiFunc( double nu0, double phi0, double * df, double nuUp0,
			double & phi, double & dphidnu0, double & dphidphi0, double & d2phidnu02, double & d2phidphi0dnu0, double & d2phidphi02 )
{

	double f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02;

	computeFourierValueAndDerivatives( nu0, phi0, df, phiNnuBasis, phiNphiBasis, nuUp0,
			f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02 );

	phi = phi0 + f;
	dphidnu0 = dfdnu0;
	dphidphi0 = 1 + dfdphi0;
	d2phidnu02 = d2fdnu02;
	d2phidphi0dnu0 = d2fdphi0dnu0;
	d2phidphi02 = d2fdphi02;

}

void FourierDeformation::muinFunc( double muin0, double nu0, double phi0, double * cf, double dmuin0dnu0, double dmuin0dphi0, double nuUp0,
			double & muin, double & dmuindnu0, double & dmuindphi0 )
{

	double f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02;

	computeFourierValueAndDerivatives( nu0, phi0, cf, muinNnuBasis, muinNphiBasis, nuUp0,
			f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02 );

	muin = muin0 + f;
	dmuindnu0 = dmuin0dnu0 + dfdnu0;
	dmuindphi0 = dmuin0dphi0 + dfdphi0;

}



double fcfuncFourier( double muin, double fc0, double acube, double R, double cnusq)
{
	double chmuin, fc;

	chmuin = cosh(muin);

	fc = fc0 - acube*R*( chmuin *(THIRD * pow(chmuin,2) - cnusq) - (THIRD - cnusq) );

	// cant contract past fc0, which is contraction down to a single line
	if ( fc > fc0 ){ fc = fc0; }

	return fc;

}



double muFuncFourier( double K, double sinRatio, double nu, double acube, double chmu0, double chmu0sq, double cnu0sq, double fc  )
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




void fcDerivativesFourier( double  muin, double dmuindnu0, double dmuindphi0,
		double dfc0dnu0, double dfc0dphi0, double acube, double R, double Rnu0, double Rphi0, double cnu,
		double cnusq, double snu, double dnudnu0, double dnudphi0, double & dfcdnu0, double & dfcdphi0 )
{

	double W, dWdnu0, dWdphi0;
	double chmuin, chmuinsq, shmuin;


	chmuin = cosh(muin);
	chmuinsq = pow(chmuin,2);
	shmuin = sinh(muin);


	// aux variable W and derivatives
	W = chmuin*( THIRD*chmuinsq - cnusq ) - (THIRD - cnusq);
	dWdnu0 = shmuin*dmuindnu0*( chmuinsq - cnusq ) + 2*cnu*snu*dnudnu0*( chmuin - 1 );
	dWdphi0 = shmuin*dmuindphi0*( chmuinsq - cnusq ) + 2*cnu*snu*dnudphi0*( chmuin - 1 );


	// final fc derivatives
	dfcdnu0 = dfc0dnu0 - acube*( Rnu0*W + R*dWdnu0 );
	dfcdphi0 = dfc0dphi0 - acube*( Rphi0*W + R*dWdphi0 );

}


/*
 *
 *
 * This is the version where the base is fixed
 *










// Main function
// first are the coordinate inputs
// second are the precomputed parameters
// third, the outputs required by the regular program
void FourierDeformation::computeDeformation( double muin0, double nu0, double phi0, double * q,
	double nuUp0, double fc0, double acube, double chmu0, double shmu0, double chmu0sq, double cnu0sq, double cphi0, double sphi0,
	double snu0, double snu0sq, double cnu0, double dfc0dnu0, double dfc0dphi0, double dmuin0dnu0, double dmuin0dphi0,
	double &mu, double &nu, double &phi, double &shmu, double &chmu, double &shmusq, double &chmusq, double &snu, double &cnu, double &snusq, double &cnusq,
	double &fc, double &dfcdnu0, double &dfcdphi0, double &K, double &sinRatio, double & dnudnu0, double &dnudphi0, double &dphidnu0, double &dphidphi0,
	double &d2nudnu02, double &d2nudphi0dnu0, double & d2nudphi02, double &d2phidnu02, double &d2phidphi02, double &d2phidphi0dnu0,
	double &dmudmu0, double &dmudnu0, double & dmudphi0 )
{



	double muin, dmuindnu0, dmuindphi0;

	// These are almost the same as c, d, and e
	// but have one extra zero for d, e
	double *coef_cf = new double[ Nc ];
	double *coef_df = new double[ Nd ];
	double *coef_ef = new double[ Ne + 1];
	pullq_parameters( q, coef_cf, coef_df, coef_ef );


	// local variables
	double R, Knu0, Kphi0, A, Bnu0, Bphi0, Rnu0, Rphi0;

	//--------------------------------------------------------------------//
	// 1. compute nu/phi and derivatives
	//--------------------------------------------------------------------//

	// Compute nu and nu derivatives
	nuFunc( nu0, phi0, coef_ef, nuUp0, nu, dnudnu0,  dnudphi0,  d2nudnu02,  d2nudphi0dnu0,  d2nudphi02 );

	// Compute phi and phi derivatives
	phiFunc( nu0, phi0, coef_df, nuUp0, phi, dphidnu0, dphidphi0, d2phidnu02, d2phidphi0dnu0, d2phidphi02 );

	//--------------------------------------------------------------------//
	// 2. compute mu
	//--------------------------------------------------------------------//

	// Need a function for this because of the singularity at nu = pi
	sinRatio = computeSinRatioFourier( nu0, nu );

	// compute K, R
	K = dnudnu0*dphidphi0 - dnudphi0*dphidnu0;
	R = sinRatio*K;

	// compute functions of nu/phi
	snu = sin(nu);
	cnu = cos(nu);
	snusq = pow(snu,2);
	cnusq = pow(cnu,2);

	// The contraction function
	muinFunc( muin0, nu0, phi0, coef_cf, dmuin0dnu0, dmuin0dphi0,  nuUp0, muin, dmuindnu0, dmuindphi0 );
	fc = fcfuncFourier( muin, fc0, acube, R, cnusq);

	// Compute mu
	mu = muFuncFourier( K, sinRatio, nu, acube, chmu0, chmu0sq, cnu0sq, fc  );


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
	fcDerivativesFourier( muin, dmuindnu0, dmuindphi0, dfc0dnu0, dfc0dphi0, acube,
			R,  Rnu0, Rphi0,cnu, cnusq, snu, dnudnu0, dnudphi0, dfcdnu0, dfcdphi0 );

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


	delete [] coef_cf;
	delete [] coef_df;
	delete [] coef_ef;

	//		cout << "nu0 = " << nu0 << ", phi0 =" << phi0 << ", mu0 = " << asinh(shmu0) << ", mu = " << mu << ", fc = " << fc << ", fc0 = " << fc0 << endl;
}

void computeFourierValueAndDerivatives( double nu0, double phi0, double * coef, int Nnu0, int Nphi0, double nuUp0,
		double & f, double & dfdnu0, double & dfdphi0, double &d2fdnu02, double &d2fdphi0dnu0, double &d2fdphi02 )
{


	double fnu0, fcphi0, fsphi0;
	double dfnu0_dnu0, dfcphi0_dphi0, dfsphi0_dphi0;
	double d2fnu0_dnu02, d2fcphi0_dphi02, d2fsphi0_dphi02;


	int k = 0;

	dfdnu0 = 0;
	dfdphi0 = 0;

	d2fdnu02 = 0;
	d2fdphi0dnu0 = 0;
	d2fdphi02 = 0;

	f = coef[k]*(1 - pow( nu0 - PI, 3 )/pow(nuUp0-PI,3));
	dfdnu0 =coef[k]*( - 3*pow( nu0 - PI, 2)/pow(nuUp0-PI,3));
	d2fdnu02 =coef[k]*( - 6*( nu0 - PI)/pow(nuUp0-PI,3));
	k++;

	// 1. Constant phi0 modes with variations in nu0
	for (int i = 0; i < Nnu0; i++)
	{
		int jcos = 0;
		int jsin = 0;
		computeFourierFunctionValuesAndDerivativesSinglePoint( nu0, phi0, i, jcos, jsin, nuUp0,
				fnu0, dfnu0_dnu0, d2fnu0_dnu02, fcphi0, dfcphi0_dphi0, d2fcphi0_dphi02,	fsphi0, dfsphi0_dphi0, d2fsphi0_dphi02 );

		f = f + coef[k]*fnu0*fcphi0;

		dfdnu0 = dfdnu0 + coef[k]*dfnu0_dnu0*fcphi0;
		dfdphi0 = dfdphi0 + coef[k]*fnu0*dfcphi0_dphi0;

		d2fdnu02 = d2fdnu02 + coef[k]*d2fnu0_dnu02*fcphi0;
		d2fdphi0dnu0 = d2fdphi0dnu0 + coef[k]*dfnu0_dnu0*dfcphi0_dphi0;
		d2fdphi02 = d2fdphi02 + coef[k]*fnu0*d2fcphi0_dphi02;

		k++;
	}


	// 2. Cycle through cos/sin modes
	for (int j = 0; j < Nphi0; j++)
	{
		for (int i = 0; i < Nnu0; i++)
		{

			int jcos = j+1;
			int jsin = j+1;
			computeFourierFunctionValuesAndDerivativesSinglePoint( nu0, phi0, i, jcos, jsin, nuUp0,
					fnu0, dfnu0_dnu0, d2fnu0_dnu02, fcphi0, dfcphi0_dphi0, d2fcphi0_dphi02,	fsphi0, dfsphi0_dphi0, d2fsphi0_dphi02 );

		    // cosine mode
			f = f + coef[k]*fnu0*fcphi0;

			dfdnu0 = dfdnu0 + coef[k]*dfnu0_dnu0*fcphi0;
			dfdphi0 = dfdphi0 + coef[k]*fnu0*dfcphi0_dphi0;

			d2fdnu02 = d2fdnu02 + coef[k]*d2fnu0_dnu02*fcphi0;
			d2fdphi0dnu0 = d2fdphi0dnu0 + coef[k]*dfnu0_dnu0*dfcphi0_dphi0;
			d2fdphi02 = d2fdphi02 + coef[k]*fnu0*d2fcphi0_dphi02;
			k++;

			// sine mode
			f = f + coef[k]*fnu0*fsphi0;

			dfdnu0 = dfdnu0 + coef[k]*dfnu0_dnu0*fsphi0;
			dfdphi0 = dfdphi0 + coef[k]*fnu0*dfsphi0_dphi0;

			d2fdnu02 = d2fdnu02 + coef[k]*d2fnu0_dnu02*fsphi0;
			d2fdphi0dnu0 = d2fdphi0dnu0 + coef[k]*dfnu0_dnu0*dfsphi0_dphi0;
			d2fdphi02 = d2fdphi02 + coef[k]*fnu0*d2fsphi0_dphi02;
			k++;


		}
	}


}


void computeFourierFunctionValuesAndDerivativesSinglePoint( double nu0, double phi0, int i, int jcos, int jsin, double nuUp0,
		double & fnu0, double &  dfnu0_dnu0, double & d2fnu0_dnu02,
		double & fcphi0, double & dfcphi0_dphi0, double & d2fcphi0_dphi02,
		double & fsphi0, double & dfsphi0_dphi0, double & d2fsphi0_dphi02 )
{
	double a = (PI*(i+1))/(PI-nuUp0);
	fnu0 = sin( a*(nu0-PI) );
	dfnu0_dnu0 = a* cos( a*(nu0-PI));
	d2fnu0_dnu02 = -a*a*sin( a*(nu0 - PI) );


	fcphi0 = cos(jcos*phi0);
	dfcphi0_dphi0 = -sin(jcos*phi0)*jcos;
	d2fcphi0_dphi02 = -cos(jcos*phi0)*jcos*jcos;

	fsphi0 = sin(jsin*phi0);
	dfsphi0_dphi0 = cos(jsin*phi0)*jsin;
	d2fsphi0_dphi02 = -sin(jsin*phi0)*jsin*jsin;

}


// Pulls the values of q
void FourierDeformation::pullq_parameters( double *q, double *coef_cf, double *coef_df, double *coef_ef )
{

	//
	for (int k = 0; k < Nc; k++ )
	{
		coef_cf[k] = q[k];
	}

	for (int k = 0; k < Nd; k++ )
		coef_df[k] = q[k+Nc];


	// set the first value to zero, and then fill in the remaining values
	// because delta nu is zero at the apex
	coef_ef[0] = 0.0;
	for (int k = 0; k < Ne; k++ )
		coef_ef[k+1] = q[k+ Nc + Nd];

}



void FourierDeformation::setSize( int muinNnuBasisIn, int muinNphiBasisIn, int nuNnuBasisIn,
		int nuNphiBasisIn, int phiNnuBasisIn, int phiNphiBasisIn )
{

    // normal version
	// MU / Nc
	muinNnuBasis = muinNnuBasisIn;
	muinNphiBasis = muinNphiBasisIn;
	Nc = muinNnuBasis*( muinNphiBasis*2+1 ) + 1;

	// NU / Ne
	nuNnuBasis = nuNnuBasisIn;
	nuNphiBasis = nuNphiBasisIn;
	Ne = nuNnuBasis*( nuNphiBasis*2+1 );

	// PHI / Nd
	phiNnuBasis = phiNnuBasisIn;
	phiNphiBasis = phiNphiBasisIn;
	Nd = phiNnuBasis*( phiNphiBasis*2+1 ) + 1;

}


void FourierDeformation::computeNParams( int & NqOut, int & NcOut, int & NdOut, int & NeOut )
{

	NcOut = Nc;
	NdOut = Nd;
	NeOut = Ne;

	NqOut = Nc+Nd+Ne;

}

FourierDeformation::FourierDeformation()
{
	setSize( 1, 1, 1, 1, 1, 1 );
}

void FourierDeformation::printNModesToConsole()
{
	cout << "Nc, Nd, Ne = " << Nc << ", " << Nd << ", " << Ne << endl;

	cout << " muinNnuBasis = " << muinNnuBasis << endl;
	cout << " muinNphiBasis = " << muinNphiBasis << endl;

	cout << " nuNnuBasis = " << nuNnuBasis << endl;
	cout << " nuNphiBasis = " << nuNphiBasis << endl;

	cout << " phiNnuBasis = " << phiNnuBasis << endl;
	cout << " phiNphiBasis = " << phiNphiBasis << endl;

}





FourierDeformation::FourierDeformation( int muinNnuBasisIn, int muinNphiBasisIn, int nuNnuBasisIn,
		int nuNphiBasisIn, int phiNnuBasisIn, int phiNphiBasisIn )
{

	setSize( muinNnuBasisIn, muinNphiBasisIn, nuNnuBasisIn,
		 nuNphiBasisIn, phiNnuBasisIn, phiNphiBasisIn );

}


void FourierDeformation::nuFunc( double nu0, double phi0, double * ef, double nuUp0,
			double & nu, double & dnudnu0, double & dnudphi0, double & d2nudnu02, double & d2nudphi0dnu0, double & d2nudphi02 )
{

	double f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02;

	computeFourierValueAndDerivatives( nu0, phi0, ef, nuNnuBasis, nuNphiBasis, nuUp0,
			f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02 );

	nu = nu0 + f;
	dnudnu0 = 1 + dfdnu0;
	dnudphi0 = dfdphi0;
	d2nudnu02 = d2fdnu02;
	d2nudphi0dnu0 = d2fdphi0dnu0;
	d2nudphi02 = d2fdphi02;

}



double computeSinRatioFourier( double nu0, double nu ){

	double sinRatio = sin(nu)/sin(nu0);

	// The sin ratio is singular at nu0 = pi, but by the definitions of
	// the transformation the ratio is 1
	if( fabs(nu0 - PI) < 1e-7)
	{
		sinRatio = 1;
	}

	return sinRatio;
}



void FourierDeformation::phiFunc( double nu0, double phi0, double * df, double nuUp0,
			double & phi, double & dphidnu0, double & dphidphi0, double & d2phidnu02, double & d2phidphi0dnu0, double & d2phidphi02 )
{

	double f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02;

	computeFourierValueAndDerivatives( nu0, phi0, df, phiNnuBasis, phiNphiBasis, nuUp0,
			f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02 );

	phi = phi0 + f;
	dphidnu0 = dfdnu0;
	dphidphi0 = 1 + dfdphi0;
	d2phidnu02 = d2fdnu02;
	d2phidphi0dnu0 = d2fdphi0dnu0;
	d2phidphi02 = d2fdphi02;

}

void FourierDeformation::muinFunc( double muin0, double nu0, double phi0, double * cf, double dmuin0dnu0, double dmuin0dphi0, double nuUp0,
			double & muin, double & dmuindnu0, double & dmuindphi0 )
{

	double f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02;

	computeFourierValueAndDerivatives( nu0, phi0, cf, muinNnuBasis, muinNphiBasis, nuUp0,
			f, dfdnu0, dfdphi0, d2fdnu02, d2fdphi0dnu0, d2fdphi02 );

	muin = muin0 + f;
	dmuindnu0 = dmuin0dnu0 + dfdnu0;
	dmuindphi0 = dmuin0dphi0 + dfdphi0;

}



double fcfuncFourier( double muin, double fc0, double acube, double R, double cnusq)
{
	double chmuin, fc;

	chmuin = cosh(muin);

	fc = fc0 - acube*R*( chmuin *(THIRD * pow(chmuin,2) - cnusq) - (THIRD - cnusq) );

	// cant contract past fc0, which is contraction down to a single line
	if ( fc > fc0 ){ fc = fc0; }

	return fc;

}



double muFuncFourier( double K, double sinRatio, double nu, double acube, double chmu0, double chmu0sq, double cnu0sq, double fc  )
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




void fcDerivativesFourier( double  muin, double dmuindnu0, double dmuindphi0,
		double dfc0dnu0, double dfc0dphi0, double acube, double R, double Rnu0, double Rphi0, double cnu,
		double cnusq, double snu, double dnudnu0, double dnudphi0, double & dfcdnu0, double & dfcdphi0 )
{

	double W, dWdnu0, dWdphi0;
	double chmuin, chmuinsq, shmuin;


	chmuin = cosh(muin);
	chmuinsq = pow(chmuin,2);
	shmuin = sinh(muin);


	// aux variable W and derivatives
	W = chmuin*( THIRD*chmuinsq - cnusq ) - (THIRD - cnusq);
	dWdnu0 = shmuin*dmuindnu0*( chmuinsq - cnusq ) + 2*cnu*snu*dnudnu0*( chmuin - 1 );
	dWdphi0 = shmuin*dmuindphi0*( chmuinsq - cnusq ) + 2*cnu*snu*dnudphi0*( chmuin - 1 );


	// final fc derivatives
	dfcdnu0 = dfc0dnu0 - acube*( Rnu0*W + R*dWdnu0 );
	dfcdphi0 = dfc0dphi0 - acube*( Rphi0*W + R*dWdphi0 );

}


*/
