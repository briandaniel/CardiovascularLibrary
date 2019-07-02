/*
 * stress.cpp
 *
 *  Created on: Sep 14, 2017
 *      Author: brian
 */

#include "stress.hpp"


void NodeGdmLV::computeStressesTractions( ParamGdmLV * prm, double t, double At ){

	elasticStress( prm );
	viscousStress( prm );

	activeStress0( prm, At );
	activeStressTD( prm, At );
	surfaceTractions( prm );
	computeTopSurfaceNorm( prm );
}


void NodeGdmLV::computeStaticStressesTractions( ParamGdmLV * prm, double At ){

	elasticStress( prm );
	activeStress0( prm, At );
	activeStressTD( prm, At );
	surfaceTractions( prm );

}

void NodeGdmLV::computeTopSurfaceNorm( ParamGdmLV * prm )
{
	double * v1vec = new double[3];
	double * v2vec = new double[3];

	double sphi = sin(phi);
	double cphi = cos(phi);

	// Compute derivatives
	double dxdmu0 = prm->a*( chmu*snu*cphi*dmudmu0 );
	double dydmu0 = prm->a*( chmu*snu*sphi*dmudmu0 );
	double dzdmu0 = prm->a*( shmu*cnu*dmudmu0 );

	double dxdphi0 = prm->a*( chmu*snu*cphi*dmudphi0 + shmu*cnu*cphi*dnudphi0 - shmu*snu*sphi*dphidphi0 );
	double dydphi0 = prm->a*( chmu*snu*sphi*dmudphi0 + shmu*cnu*sphi*dnudphi0 + shmu*snu*cphi*dphidphi0 );
	double dzdphi0 = prm->a*( shmu*cnu*dmudphi0 - chmu*snu*dnudphi0 );

	// Assign vectors
	v1vec[0] = dxdmu0;
	v1vec[1] = dydmu0;
	v1vec[2] = dzdmu0;

	v2vec[0] = dxdphi0;
	v2vec[1] = dydphi0;
	v2vec[2] = dzdphi0;


	// Compute cross product
	crossProduct3( v1vec, v2vec, vcrossTop);
	vCrossTopNorm = vectorNorm (vcrossTop, 3, 2 );

	// unit normal
	for(int i = 0; i < 3; i++)
		nTop[i] = vcrossTop[i]/vCrossTopNorm;

	delete [] v1vec;
	delete [] v2vec;


}

void NodeGdmLV::surfaceTractions( ParamGdmLV * prm )
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
	crossProduct3( v1, v2, vcrossSide );

	// norm of the side vectors
	vCrossSideNorm = vectorNorm(vcrossSide,3,2);

	// unit normal
	for(int i = 0; i < 3; i++)
		nSide[i] = vcrossSide[i]/vCrossSideNorm;



	// clean up
	delete [] v1;
	delete [] v2;


}




void NodeGdmLV::elasticStress( ParamGdmLV * prm ){




	// temporary variables
	double ** Sefib = new double * [3];
	for(int i = 0; i < 3; i ++){ Sefib[i] = new double [3]; }

	double eW;

	eW = exp(   prm->bff*pow(Efib[2][2],2)
	          + prm->bxx*( pow(Efib[1][1],2) + pow(Efib[0][0],2) + pow(Efib[0][1],2) + pow(Efib[1][0],2) )
			  + prm->bfx*( pow(Efib[0][2],2) + pow(Efib[2][0],2) + pow(Efib[1][2],2) + pow(Efib[2][1],2) ) );

	Sefib[0][0] = prm->ke*eW*prm->bxx*Efib[0][0];
	Sefib[0][1] = prm->ke*eW*prm->bxx*Efib[0][1];
	Sefib[1][1] = prm->ke*eW*prm->bxx*Efib[1][1];

	Sefib[0][2] = prm->ke*eW*prm->bfx*Efib[0][2];
	Sefib[1][2] = prm->ke*eW*prm->bfx*Efib[1][2];

	Sefib[2][2] = prm->ke*eW*prm->bff*Efib[2][2];

	// symmetric
	Sefib[1][0] = Sefib[0][1];
	Sefib[2][0] = Sefib[0][2];
	Sefib[2][1] = Sefib[1][2];

	rotateFiberToProlate( Sefib, Se );


	for(int i = 0; i < 3; i ++){ delete [] Sefib[i]; }
	delete [] Sefib;

	// only use this if set to use a rivlin-mooney constitutive law
	// in this case bff = c1 and bxx = c2. bfx, ke are unused.
	if(prm->rivlinMooneyMaterial == 1)
	{
		elasticStressRivlinMooney( prm );
	}


}



void NodeGdmLV::elasticStressRivlinMooney( ParamGdmLV * prm )
{

	Se[0][0] = 2*( prm->bff + 2*prm->bxx*( E[1][1] + E[2][2] + 1 ) ) ;
	Se[1][1] = 2*( prm->bff + 2*prm->bxx*( E[0][0] + E[2][2] + 1 ) ) ;
	Se[2][2] = 2*( prm->bff + 2*prm->bxx*( E[0][0] + E[1][1] + 1 ) ) ;

	Se[0][1] = -4*prm->bxx*E[0][1];
	Se[0][2] = -4*prm->bxx*E[0][2];
	Se[1][2] = -4*prm->bxx*E[1][2];

	Se[1][0] = Se[0][1];
	Se[2][0] = Se[0][2];
	Se[2][1] = Se[1][2];
}


void NodeGdmLV::viscousStress( ParamGdmLV * prm ){

	// temp matrix
	double ** Atemp = new double * [3];
	for (int k = 0; k < 3; k++)
	{
		Atemp[k] = new double [3];
	}


	for( int k = 0; k < prm->Nq; k++ ){

		// Svqk is computed using matrix multiplication
		matrixMultiply( Cinv, dEdq[k], Atemp, 3, 3, 3);
		matrixMultiply( Atemp, Cinv, Svq[k], 3, 3, 3);

		// There is also a constant factor; multiply by it
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				Svq[k][i][j] = 2*prm->kv*Svq[k][i][j];
			}
		}

	}


	// clean up
	for (int k = 0; k < 3; k++)
	{
		delete [] Atemp[k];
	}
	delete [] Atemp;
}

// active stress time derivative componenent
void NodeGdmLV::activeStress0( ParamGdmLV * prm, double At )
{

	// local definitions
	double G, Tff;
	double ** Sa0fib = new double * [3];
	for(int k = 0; k < 3; k ++){ Sa0fib[k] = new double [3]; }


	// compute stretch ratio in fiber direction
	lambda = sqrt( 2*Efib[2][2] + 1 );

	// Length-tension factor
	G = tensionFunc( prm, lambda );

	// Fiber direction tension
	Tff = At*prm->ka*G;

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


	// clean up
	for(int k = 0; k < 3; k ++){ delete [] Sa0fib[k]; }
	delete [] Sa0fib;


}



// active stress component not multiplied by a time derivative
void NodeGdmLV::activeStressTD( ParamGdmLV * prm, double At )
{



	// local definitions
	double G;

	double * Tfibq = new double [prm->Nq];

	double ** Safibqi = new double * [3];
	for(int k = 0; k < 3; k ++){ Safibqi[k] = new double [3]; }


	// compute stretch ratio in fiber direction
	lambda = sqrt( 2*Efib[2][2] + 1 );

	// Length-tension factor
	G = tensionFunc( prm, lambda );


	// q dependent terms
	for (int k = 0; k < prm->Nq; k++)
	{
		// Compute the derivative of lambda
		double dlambda_dq = (1.0/lambda)*dEfibdq[k][2][2];

		// Fiber direction tension derivatives
		Tfibq[k] = At*G*prm->kav*dlambda_dq;
	}

	// Compute the stress tensor qi derivatives
	for (int k = 0; k < prm->Nq; k++)
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



	// clean up
	for(int k = 0; k < 3; k ++){ delete [] Safibqi[k]; }
	delete [] Safibqi;
	delete [] Tfibq;



}



double tensionFunc( ParamGdmLV * prm, double lambda)
{
	double LS;
	double geff;

	LS = prm->Ls0*lambda;

	geff = exp( - pow( LS - prm->Lsmax, 2) / ( prm->twoLswSq ) );

	return geff;

}

double activationFunc( double t, double eff_ed, ParamGdmLV * prm )
{
	double At = 0;
	double t_cycle = mod(t, prm->Tc);
	double d = 1/(1+prm->kd*eff_ed);


	if( t_cycle < (prm->Tc-prm->Ta) )
	{
		At = 0;
	}else{
		At = pow( sin(PI* (t_cycle - prm->Tc + prm->Ta )/prm->Ta ), d*prm->activePower);
	}


	return At;


}



// assumes that F, FT, Se, Sa0 have been computed
void NodeGdmLV::computeCauchyStress( )
{

	// temp matrix
	double ** Stemp = new double * [3];
	double ** Atemp = new double * [3];
	double ** QprolateToCart = new double * [3];
	double ** QprolateToCartT = new double * [3];

	for (int k = 0; k < 3; k++)
	{
		Atemp[k] = new double [3];
		Stemp[k] = new double [3];
		QprolateToCart[k] = new double [3];
		QprolateToCartT[k] = new double [3];

	}


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
	matrixMultiply( F, Stemp, Atemp, 3, 3, 3);
	matrixMultiply( Atemp, FT, sigma_prol, 3, 3, 3);



	// Compute cauchy stress in cartesian coordinates
	computeProlateToCartesianRotationMatrix( mu, nu, phi, QprolateToCart );
	matrixTranspose( QprolateToCart, QprolateToCartT, 3, 3);

	matrixMultiply( QprolateToCart, sigma_prol, Atemp, 3, 3, 3 );
	matrixMultiply( Atemp, QprolateToCartT, sigma_cart, 3, 3, 3 );


	// clean up
	for (int k = 0; k < 3; k++)
	{
		delete [] Stemp[k];
		delete [] Atemp[k];
		delete [] QprolateToCart[k];
		delete [] QprolateToCartT[k];


	}
	delete [] Stemp;
	delete [] Atemp;
	delete [] QprolateToCart;
	delete [] QprolateToCartT;

}

























