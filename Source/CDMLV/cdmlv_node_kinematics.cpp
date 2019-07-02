/*
 * cdmlv_node_kinematics.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */

#include "cdmlv_node.hpp"



void CDMNode::computeDeformation( vector <double> & q, ParamCDMLV & prm )
{
	// Right now this is the only deformation model
	computeFourierDeformation( q, prm );

	// Compute the cartesian coordinates, etc.
	computeSubvalues();

}



void CDMNode::computeFourierDeformation( vector <double> & q, ParamCDMLV & prm )
{

	// use the LV deformation function library to compute the full deformation
	prm.fourierDef.computeDeformation( muin0, nu0, phi0, q.data(), nuUp0, fc0, acube, chmu0, shmu0, chmu0sq, cnu0sq, cphi0, sphi0,
		snu0, snu0sq, cnu0, dfc0dnu0, dfc0dphi0, dmuin0dnu0, dmuin0dphi0,
		mu, nu, phi, shmu, chmu, shmusq, chmusq, snu, cnu, snusq, cnusq,
		fc, dfcdnu0, dfcdphi0, K, sinRatio,  dnudnu0, dnudphi0, dphidnu0, dphidphi0,
		d2nudnu02, d2nudphi0dnu0,  d2nudphi02, d2phidnu02, d2phidphi02, d2phidphi0dnu0,
		dmudmu0, dmudnu0,  dmudphi0 );

}


void CDMNode::computeSubvalues( )
{

	prolate2xyz( mu, nu, phi, a, x, y, z);

	gmu = a*sqrt( shmusq + snusq );
	gnu = gmu;
	gphi = a*shmu*snu;

	// only need the displacement on the surface nodes
	udisp[0] = x - x0;
	udisp[1] = y - y0;
	udisp[2] = z - z0;

}



void CDMNode::computeDeformationGradient()
{

	F[0][0] = gmu/gmu0*dmudmu0;
	F[1][0] = 0;
	F[2][0] = 0;

	F[0][1] = ( gmu/gnu0 )*dmudnu0;
	F[1][1] = ( gnu/gnu0 )*dnudnu0;
	F[2][1] = ( gphi/gnu0 )*dphidnu0;

	F[0][2] = ( gmu/gphi0 )*dmudphi0;
	F[1][2] = ( gnu/gphi0 )*dnudphi0;
	F[2][2] = ( gphi/gphi0 )*dphidphi0;


	// At the apex the values must be adjusted
	if( fabs(nu0 - PI) < 1e-7)
	{

		F[0][2] = 0;
		F[1][2] = 0;
		F[2][2] = ( shmu / shmu0 )*( dphidphi0 );
	}

	computeMatrixInverse33(F, Finv);
	matrixTranspose(F, FT );

}


void CDMNode::computeStrainTensors()
{
	vector< vector<double> > FTinv (3, vector<double>(3,0));
	vector< vector<double> > Binv (3, vector<double>(3,0));

	// compute inverse of F^T
	computeMatrixInverse33( FT, FTinv );

	// compute C, Binv
	matrixMultiply( FT, F, C );
	matrixMultiply( FTinv, Finv, Binv );
	matrixMultiply( Finv, FTinv, Cinv );

	// Compute E
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{

			if( i == j)
			{
				E[i][j] = 0.5*(C[i][j] - 1);
			}else{
				E[i][j] = 0.5*C[i][j];
			}
		}
	}

	// Compute e
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{

			if( i == j)
			{
				estrain[i][j] = 0.5*( 1 - Binv[i][j] );
			}else{
				estrain[i][j] = -0.5*Binv[i][j];
			}
		}
	}

	rotateProlateToFiber( E, Efib );
	rotateProlateToFiber( estrain, efib );



}



void CDMNode::rotateProlateToFiber( vector< vector<double> > & T, vector< vector<double> > & Tfib )
{

	vector< vector<double> > Atemp (3, vector<double>(3,0));

	// multiply to compute rotation
	matrixMultiply( QT, T, Atemp );
	matrixMultiply( Atemp, Q, Tfib );

}



void CDMNode::rotateFiberToProlate( vector< vector<double> > & Tfib, vector< vector<double> > &  T )
{

	vector< vector<double> > Atemp (3, vector<double>(3,0));

	// multiply to compute rotation
	matrixMultiply( Q, Tfib, Atemp );
	matrixMultiply( Atemp, QT, T );

}



void CDMNode::computeCavityVolume( ParamCDMLV & prm )
{
	// This value is only used if the node is on the endocardial surface
	if( endoSurfNode == 1 )
	{
		// because this value is only computed at the endocardial surface,
		// chmu is automatically chmuin0
		VI = acube*K*snu*( 1.0/3.0*(pow(chmu,3) - 1) + pow( cnu, 2 )*(1 - chmu));

	}else{
		VI = 0;
	}


}


// compute the deformations/strains in the domain
void CDMNode::kinematicComputations( vector <double> & q, ParamCDMLV & prm )
{

		// compute deformation, strains, and positions (at the boundary) at qdqk
		computeDeformation( q, prm) ;

		computeDeformationGradient();

		computeStrainTensors();

		computeCavityVolume( prm );

}








