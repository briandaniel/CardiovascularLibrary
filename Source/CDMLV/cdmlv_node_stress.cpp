/*
 * cdmlv_node_stress.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */

#include "cdmlv_node.hpp"




void CDMNode::computeStressesTractions( ParamCDMLV & prm, double t, double At ){

	elasticStress( prm );
	viscousStress( prm );

	activeStress0( prm, At );
	activeStressTD( prm, At );
	surfaceTractions( prm );
	computeTopSurfaceNorm( prm );
}


void CDMNode::computeStaticStressesTractions( ParamCDMLV & prm, double At ){

	elasticStress( prm );
	activeStress0( prm, At );
	activeStressTD( prm, At );
	surfaceTractions( prm );

}

void CDMNode::computeTopSurfaceNorm( ParamCDMLV & prm )
{
	double * v1vec = new double[3];
	double * v2vec = new double[3];

	double sphi = sin(phi);
	double cphi = cos(phi);

	// Compute derivatives
	double dxdmu0 = prm.a*( chmu*snu*cphi*dmudmu0 );
	double dydmu0 = prm.a*( chmu*snu*sphi*dmudmu0 );
	double dzdmu0 = prm.a*( shmu*cnu*dmudmu0 );

	double dxdphi0 = prm.a*( chmu*snu*cphi*dmudphi0 + shmu*cnu*cphi*dnudphi0 - shmu*snu*sphi*dphidphi0 );
	double dydphi0 = prm.a*( chmu*snu*sphi*dmudphi0 + shmu*cnu*sphi*dnudphi0 + shmu*snu*cphi*dphidphi0 );
	double dzdphi0 = prm.a*( shmu*cnu*dmudphi0 - chmu*snu*dnudphi0 );

	// Assign vectors
	v1vec[0] = dxdmu0;
	v1vec[1] = dydmu0;
	v1vec[2] = dzdmu0;

	v2vec[0] = dxdphi0;
	v2vec[1] = dydphi0;
	v2vec[2] = dzdphi0;


	// Compute cross product
	crossProduct3( v1vec, v2vec, vcrossTop.data());
	vCrossTopNorm = vectorNorm (vcrossTop.data(), 3, 2 );

	// unit normal
	for(int i = 0; i < 3; i++)
		nTop[i] = vcrossTop[i]/vCrossTopNorm;

	delete [] v1vec;
	delete [] v2vec;


}

void CDMNode::surfaceTractions( ParamCDMLV & prm )
{

	double cphi, sphi;

	double * v1 = new double [3];
	double * v2 = new double [3];

	sphi = sin(phi);
	cphi = cos(phi);

	// Compute derivatives
	dxdnu0 = a*( chmu*snu*cphi*dmudnu0 + shmu*cnu*cphi*dnudnu0 - shmu*snu*sphi*dphidnu0 );
	dydnu0 = a*( chmu*snu*sphi*dmudnu0 + shmu*cnu*sphi*dnudnu0 + shmu*snu*cphi*dphidnu0 );
	dzdnu0 = a*( shmu*cnu*dmudnu0 - chmu*snu*dnudnu0 );

	dxdphi0 = a*( chmu*snu*cphi*dmudphi0 + shmu*cnu*cphi*dnudphi0 - shmu*snu*sphi*dphidphi0 );
	dydphi0 = a*( chmu*snu*sphi*dmudphi0 + shmu*cnu*sphi*dnudphi0 + shmu*snu*cphi*dphidphi0 );
	dzdphi0 = a*( shmu*cnu*dmudphi0 - chmu*snu*dnudphi0 );

	// Assign vectors
	v1[0] = dxdnu0;
	v1[1] = dydnu0;
	v1[2] = dzdnu0;

	v2[0] = dxdphi0;
	v2[1] = dydphi0;
	v2[2] = dzdphi0;

	// Compute cross product
	crossProduct3( v1, v2, vcrossSide.data() );

	// norm of the side vectors
	vCrossSideNorm = vectorNorm(vcrossSide.data(),3,2);

	// unit normal
	for(int i = 0; i < 3; i++)
		nSide[i] = vcrossSide[i]/vCrossSideNorm;



	// clean up
	delete [] v1;
	delete [] v2;


}




void CDMNode::elasticStress( ParamCDMLV & prm ){

	// temporary variables
	vector< vector<double> > Sefib (3, vector<double>(3,0));
	double eW;

	eW = exp(   prm.bff*pow(Efib[2][2],2)
	          + prm.bxx*( pow(Efib[1][1],2) + pow(Efib[0][0],2) + pow(Efib[0][1],2) + pow(Efib[1][0],2) )
			  + prm.bfx*( pow(Efib[0][2],2) + pow(Efib[2][0],2) + pow(Efib[1][2],2) + pow(Efib[2][1],2) ) );

	Sefib[0][0] = prm.ke*eW*prm.bxx*Efib[0][0];
	Sefib[0][1] = prm.ke*eW*prm.bxx*Efib[0][1];
	Sefib[1][1] = prm.ke*eW*prm.bxx*Efib[1][1];

	Sefib[0][2] = prm.ke*eW*prm.bfx*Efib[0][2];
	Sefib[1][2] = prm.ke*eW*prm.bfx*Efib[1][2];

	Sefib[2][2] = prm.ke*eW*prm.bff*Efib[2][2];

	// symmetric
	Sefib[1][0] = Sefib[0][1];
	Sefib[2][0] = Sefib[0][2];
	Sefib[2][1] = Sefib[1][2];

	rotateFiberToProlate( Sefib, Se );

}



void CDMNode::viscousStress( ParamCDMLV & prm ){

	// temp matrix
	vector< vector<double> > Atemp (3, vector<double>(3,0));


	for( int k = 0; k < prm.Nq; k++ ){

		// Svqk is computed using matrix multiplication
		matrixMultiply( Cinv, dEdq[k], Atemp );
		matrixMultiply( Atemp, Cinv, Svq[k] );

		// There is also a constant factor; multiply by it
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				Svq[k][i][j] = 2*prm.kv*Svq[k][i][j];
			}
		}

	}

}

// active stress time derivative componenent
void CDMNode::activeStress0( ParamCDMLV & prm, double At )
{

	// local definitions
	double G, Tff;
	vector< vector<double> > Sa0fib  (3, vector<double>(3,0));


	// compute stretch ratio in fiber direction
	lambda = sqrt( 2*Efib[2][2] + 1 );

	// Length-tension factor
	G = gaussTensionFunc( lambda, prm.Ls0, prm.Lsmax, prm.Lsw );


	// Fiber direction tension
	Tff = At*prm.ka*G;

	// Set fiber direction to be the tension and all others zero
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			Sa0fib[i][j] = 0.0;
		}
	}
	Sa0fib[2][2] = Tff;


	// Rotate back into prolate coordinates
	rotateFiberToProlate( Sa0fib, Sa0 );



}



// active stress component not multiplied by a time derivative
void CDMNode::activeStressTD( ParamCDMLV & prm, double At )
{

	// local definitions
	double G;

	vector<double> Tfibq( Nq, 0);
	vector< vector<double> > Safibqi  (3, vector<double>(3,0));


	// compute stretch ratio in fiber direction
	lambda = sqrt( 2*Efib[2][2] + 1 );

	// Length-tension factor
	G = gaussTensionFunc( lambda, prm.Ls0, prm.Lsmax, prm.Lsw );


	// q dependent terms
	for (int k = 0; k < prm.Nq; k++)
	{
		// Compute the derivative of lambda
		double dlambda_dq = (1.0/lambda)*dEfibdq[k][2][2];

		// Fiber direction tension derivatives
		Tfibq[k] = At*G*prm.kav*dlambda_dq;
	}

	// Compute the stress tensor qi derivatives
	for (int k = 0; k < prm.Nq; k++)
	{
		// Set fiber direction to be the tension and all others zero
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				Safibqi[i][j] = 0.0;
			}
		}
		Safibqi[2][2] = Tfibq[k]; // use the correct q derivative indexed by k

		// Rotate back into prolate coordinates
		// store in the correct q index: k
		rotateFiberToProlate( Safibqi, Saq[k] );

	}


}




// assumes that F, FT, Se, Sa0 have been computed
void CDMNode::computeCauchyStress()
{

	// temp matrix
	vector< vector<double> > Stemp (3, vector<double>(3,0));
	vector< vector<double> > Atemp (3, vector<double>(3,0));
	vector< vector<double> > QprolateToCart (3, vector<double>(3,0));
	vector< vector<double> > QprolateToCartT (3, vector<double>(3,0));

	// set S to the required stress components
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			Stemp[i][j] = Se[i][j]+ Sa0[i][j];
		}
	}


	// Compute cauchy stress in prolate coordinates
	// sigma_prol = F_prol S_prol F_prol^T
	matrixMultiply( F, Stemp, Atemp );
	matrixMultiply( Atemp, FT, sigma_prol );

	// Compute cauchy stress in cartesian coordinates
	computeProlateToCartesianRotationMatrix( mu, nu, phi, QprolateToCart );
	matrixTranspose( QprolateToCart, QprolateToCartT );

	matrixMultiply( QprolateToCart, sigma_prol, Atemp );
	matrixMultiply( Atemp, QprolateToCartT, sigma_cart );


}




