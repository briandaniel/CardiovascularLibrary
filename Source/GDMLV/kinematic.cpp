/*
 * kinematic.cpp
 *
 *  Created on: Sep 12, 2017
 *      Author: brian
 */


#include "kinematic.hpp"

void NodeGdmLV::computeDeformation( double * q, FourierDeformation & fourierDef, ParamGdmLV * prm )
{

	if( prm->deformationModelIndex == 1 )
	{
		computeFourierDeformation( q, fourierDef );
	}
	else
	{
		computeNodeSevenParamDeformation( q );
	}

	/*
	cout << "--------" << endl;
	cout << mu0 << " " << nu0 << " " << phi0 << endl;
	cout << mu << " " << nu << " " << phi << endl;
	*/
}


void NodeGdmLV::computeFourierDeformation( double * q, FourierDeformation & fourierDef )
{

	// use the LV deformation function library to compute the full deformation
	fourierDef.computeDeformation( muin0, nu0, phi0, q, nuUp0, fc0, acube, chmu0, shmu0, chmu0sq, cnu0sq, cphi0, sphi0,
		snu0, snu0sq, cnu0, dfc0dnu0, dfc0dphi0, dmuin0dnu0, dmuin0dphi0,
		mu, nu, phi, shmu, chmu, shmusq, chmusq, snu, cnu, snusq, cnusq,
		fc, dfcdnu0, dfcdphi0, K, sinRatio,  dnudnu0, dnudphi0, dphidnu0, dphidphi0,
		d2nudnu02, d2nudphi0dnu0,  d2nudphi02, d2phidnu02, d2phidphi02, d2phidphi0dnu0,
		dmudmu0, dmudnu0,  dmudphi0 );

}


void NodeGdmLV::computeNodeSevenParamDeformation( double * q )
{
	computeSevenParamDeformation( muin0, nu0, phi0, q, nuUp0, fc0, acube, chmu0, shmu0, chmu0sq, cnu0sq, cphi0, sphi0,
			snu0, snu0sq, cnu0, dfc0dnu0, dfc0dphi0, dmuin0dnu0, dmuin0dphi0, chmue,
			mu, nu, phi, shmu, chmu, shmusq, chmusq, snu, cnu, snusq, cnusq,
			fc, dfcdnu0, dfcdphi0, K, sinRatio,  dnudnu0, dnudphi0, dphidnu0, dphidphi0,
			d2nudnu02, d2nudphi0dnu0,  d2nudphi02, d2phidnu02, d2phidphi02, d2phidphi0dnu0,
			dmudmu0, dmudnu0,  dmudphi0 );

}



void NodeGdmLV::computeDeformationGradient()
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
	matrixTranspose(F, FT, 3, 3);

}


void NodeGdmLV::computeStrainTensors()
{
	double ** FTinv = new double * [3];
	double ** Binv = new double *[3];
	for (int k = 0; k < 3; k++)
	{
		FTinv[k] = new double [3];
		Binv[k] = new double [3];

	}

	// compute inverse of F^T
	computeMatrixInverse33( FT, FTinv );

	// compute C, Binv
	matrixMultiply( FT, F, C, 3, 3, 3);
	matrixMultiply( FTinv, Finv, Binv, 3, 3, 3);
	matrixMultiply( Finv, FTinv, Cinv, 3, 3, 3);

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

	// clean up
	for (int k = 0; k < 3; k++)
	{
		delete [] FTinv[k];
		delete [] Binv[k];
	}
	delete [] FTinv;
	delete [] Binv;

}



void NodeGdmLV::rotateProlateToFiber( double** T, double** Tfib )
{

	double ** Atemp = new double * [3];
	for (int k = 0; k < 3; k++)
	{
		Atemp[k] = new double [3];
	}

	// multiply to compute rotation
	matrixMultiply( QT, T, Atemp, 3, 3, 3);
	matrixMultiply( Atemp, Q, Tfib, 3, 3, 3);

	// clean up
	for (int k = 0; k < 3; k++)
	{
		delete [] Atemp[k];
	}
	delete [] Atemp;

}



void NodeGdmLV::rotateFiberToProlate( double** Tfib, double** T )
{

	double ** Atemp = new double * [3];
	for (int k = 0; k < 3; k++)
	{
		Atemp[k] = new double [3];
	}

	// multiply to compute rotation
	matrixMultiply( Q, Tfib, Atemp, 3, 3, 3);
	matrixMultiply( Atemp, QT, T, 3, 3, 3);

	// clean up
	for (int k = 0; k < 3; k++)
	{
		delete [] Atemp[k];
	}
	delete [] Atemp;

}




// Compute the displacements assuming that x has been computed as desired. x0, y0, z0 are fixed.
void NodeGdmLV::computeDisplacements( ){

	// only need the displacement on the surface nodes
	udisp[0] = x - x0;
	udisp[1] = y - y0;
	udisp[2] = z - z0;

}





// compute the deformations/strains in the domain
void NodeGdmLV::deformationStrainComputations( double * q, FourierDeformation & fourierDef, ParamGdmLV * prm )
{

		// compute deformation, strains, and positions (at the boundary) at qdqk
		computeDeformation( q, fourierDef, prm) ;

		computeSubvalues();

		computeDeformationGradient();

		computeStrainTensors();

		computeDisplacements();

		computeCavityVolume( prm );

}





void NodeGdmLV::computeCavityVolume( ParamGdmLV * prm )
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






