/*
 * cdmlv_lid.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */


#include "cdmlv_lid.hpp"

void LidCDMLV::setupLid( ParamCDMLV & prm )
{

	// scale of prolate coordinate system
	a = prm.a;
	nuUp0Min = prm.nuUp0Min;

	// number of deformation parameters
	Nq = prm.Nq;

	// lid size
	NmuLid = prm.NmuLid;
	Nr = NmuLid;
	Nphi = prm.Nphi;


	// vectors that give the values at the endocardial boundary at the base
	mu0Vec.assign(Nphi,0);
	nu0Vec.assign(Nphi,0);
	phi0Vec.assign(Nphi,0);

	muVec.assign(Nphi,0);
	nuVec.assign(Nphi,0);
	phiVec.assign(Nphi,0);

	nuUp0Vec.assign(Nphi,0);
	dmuin0dnu0Vec.assign(Nphi,0);
	dmuin0dphi0Vec.assign(Nphi,0);

	fc0Vec.assign(Nphi,0);
	dfc0dnu0Vec.assign(Nphi,0);
	dfc0dphi0Vec.assign(Nphi,0);
	dmudphi0Vec.assign(Nphi,0);
	dnudphi0Vec.assign(Nphi,0);
	dphidphi0Vec.assign(Nphi,0);

	// second surface vectors at the boundary
	rVec.assign(Nphi,0);
	zVec.assign(Nphi,0);


	// 2D arrays for surface 1
	dS0.assign(NmuLid,vector<double>(Nphi,0));
	mu0Surf1.assign(NmuLid,vector<double>(Nphi,0));
	nu0Surf1.assign(NmuLid,vector<double>(Nphi,0));
	phi0.assign(NmuLid,vector<double>(Nphi,0));
	muSurf1.assign(NmuLid,vector<double>(Nphi,0));
	nuSurf1.assign(NmuLid,vector<double>(Nphi,0));
	phi.assign(NmuLid,vector<double>(Nphi,0));

	zSurf1.assign(NmuLid,vector<double>(Nphi,0));
	r.assign(NmuLid,vector<double>(Nphi,0));
	zLid.assign(NmuLid,vector<double>(Nphi,0));
	ur.assign(NmuLid,vector<double>(Nphi,0));
	x0.assign(NmuLid,vector<double>(Nphi,0));
	y0.assign(NmuLid,vector<double>(Nphi,0));
	z0.assign(NmuLid,vector<double>(Nphi,0));
	x.assign(NmuLid,vector<double>(Nphi,0));
	y.assign(NmuLid,vector<double>(Nphi,0));
	z.assign(NmuLid,vector<double>(Nphi,0));

	uLid.assign(NmuLid,vector<vector<double>>(Nphi,vector<double>(3,0)));
	duLiddq.assign(Nq, vector<vector<vector<double>>>(NmuLid,vector<vector<double>>(Nphi, vector<double>(3,0))));
	vCross.assign(NmuLid,vector<vector<double>>(Nphi,vector<double>(3,0)));
	vur.assign(NmuLid,vector<vector<double>>(Nphi,vector<double>(3,0)));
	vphi0.assign(NmuLid,vector<vector<double>>(Nphi,vector<double>(3,0)));

	cylWeight.assign(NmuLid,vector<double>(Nphi,0));

	// volume derivatives
	dVirreg1dq.assign(Nq,0);
	dVirreg2dq.assign(Nq,0);

	setInitialPositions( prm );

}


void LidCDMLV::setInitialPositions( ParamCDMLV & prm ){


	acube = pow(a,3);

	// Compute phi0
    dphi0 = 2 * PI / Nphi;
    for(int k = 0; k < Nphi; k++)
    {
    	phi0Vec[k] = dphi0*k;
    	nuUp0Vec[k] = computeNuUp0( phi0Vec[k], prm.nuUp0Spline );
    	nu0Vec[k] = nuUp0Vec[k];

    	mu0Vec[k] = prm.endo0.evaluateSplines(nu0Vec[k],phi0Vec[k]);

    }



	double cmu, cphi;
	for(int i = 0; i < NmuLid; i++)
	{
		for(int j = 0; j < Nphi; j++)
		{
			cmu = computeIntegralCoef ( i, NmuLid );
			cphi = computePeriodicIntegralCoef ( j, Nphi );

			dS0[i][j] = cmu*cphi;
		}
	}


	// Compute mu0, nu0, phi0
	double nuAxis = nuUp0Min;

	for(int i = 0; i < NmuLid; i++)
	{
		for(int j = 0; j < Nphi; j++)
		{
			double dmu0 = ( mu0Vec[j] - 0.0 )/(NmuLid - 1.0);

			mu0Surf1[i][j] = dmu0*i;
			phi0[i][j] = phi0Vec[j];
			nu0Surf1[i][j] = nuAxis + ( nu0Vec[j] - nuAxis )* mu0Surf1[i][j]/mu0Vec[j];

		}
	}




	for(int k = 0; k < Nphi; k++)
	{
		double snu0 = sin(nu0Vec[k]);
		double cnu0 = cos(nu0Vec[k]);

		double snu0sq = pow(snu0,2);
		double cnu0sq = pow(cnu0,2);

		double shmu0 = sinh(mu0Vec[k]);
		double chmu0 = cosh(mu0Vec[k]);

		double shmu0sq = pow(shmu0,2);
		double chmu0sq = pow(chmu0,2);

		double chmuin0 = cosh(mu0Vec[k]);
		double chmuin0sq = pow(chmuin0,2);
		double shmuin0 = sinh(mu0Vec[k]);

		double cphi0 = cos(phi0Vec[k]);
		double sphi0 = sin(phi0Vec[k]);

		// evaluate spline surface derivatives used in additional computations (derivative of the constant part is 0)
		double dmuout0dnu0, dmuout0dphi0;
		prm.endo0.evaluateSplineDerivatives( nu0Vec[k], phi0Vec[k], dmuin0dnu0Vec[k], dmuin0dphi0Vec[k]);
		prm.epi0.evaluateSplineDerivatives( nu0Vec[k], phi0Vec[k], dmuout0dnu0, dmuout0dphi0);


		// compute the fc0 function and its derivatives
		double dfc0dnu0, dfc0dphi0;
		fc0func( acube, chmuin0, chmuin0sq, shmuin0, dmuin0dnu0Vec[k], dmuin0dphi0Vec[k],
				 cnu0, cnu0sq, snu0, cphi0, fc0Vec[k], dfc0dnu0Vec[k], dfc0dphi0Vec[k]);


	}

	// compute the initial x, y, z positions
	// dummy used to compute unstressed locations
	FourierDeformation fourierDefTemp( prm.muinNnuGrid,  prm.muinNphiGrid,  prm.nuNnuGrid,  prm.nuNphiGrid,	 prm.phiNnuGrid,  prm.phiNphiGrid );

	// compute zero pressure position of lid
	vector<double> qTemp(Nq,0);
	computeFirstSurfacePosition( qTemp, prm );

	computeSecondSurfacePosition();

	for(int i = 0; i < NmuLid; i++)
	{
		for(int j = 0; j < Nphi; j++)
		{
			x0[i][j] = x[i][j];
			y0[i][j] = y[i][j];
			z0[i][j] = z[i][j];
		}
	}

}





void LidCDMLV::computeLid( vector<double> & q, ParamCDMLV & prm )
{

	// compute the positions
	computeFirstSurfacePosition( q, prm );
	computeSecondSurfacePosition();

	// compute the volumes
	computeFirstIrregularVolume();
	computeSecondIrregularVolume();

	// compute the displacement/displacement derivatives
	computeDisplacement();

	// compute the normal surface vectors
	computeSurfaceNormals();

}

void LidCDMLV::exportVolumes( double & Virreg, vector<double> & dVirregdq )
{

	Virreg = Virreg1 + Virreg2;

	for( int m = 0; m < Nq; m++ )
	{
		dVirregdq[m] = dVirreg1dq[m] + dVirreg2dq[m];
	}

}

// assumes the lid has been computed
void LidCDMLV::computeLidSurfaceIntegrals( vector<double> & etal )
{



	double * integrand_r = new double [Nr];
	double * Ir = new double [Nphi];

	// different value of eta for each q
	for(int m = 0; m < Nq; m++)
	{


		// compute the integral in the r direction first
		// this integral has uneven spacing so it is computed
		// using the trapezoid rule with uneven spacing
		for(int j = 0; j < Nphi; j++)
		{

			// compute the integrand
			for(int i = 0; i < Nr; i++)
			{
				integrand_r[i] = dotProd( duLiddq[m][i][j], vCross[i][j]  );

			}

			// compute the integral
			Ir[j] = 0;
			for(int i = 0; i < Nr-1; i++)
			{
				double dur = ur[i+1][j] - ur[i][j];
				Ir[j] = Ir[j] + dur/2.0 * ( integrand_r[i+1] + integrand_r[i] );
			}

		}


		// compute the integral in the phi direction
		etal[m] = 0;
		for(int j = 0; j < Nphi; j++)
		{
			double cphi = computePeriodicIntegralCoef ( j, Nphi );
			etal[m] = etal[m] + cphi*Ir[j]*dphi0;
		}

	}


	// clean up
	delete [] integrand_r;
	delete [] Ir;


}


void LidCDMLV::computeSurfaceNormals()
{

	// cycle through phi first so that
	// boundary computations are computed
	// in the correct order
	for(int j = 0; j < Nphi; j++)
	{
		double drOutdphi0 = a*sin(nuVec[j])*cosh(muVec[j])*dmudphi0Vec[j]
						    + a*sinh(muVec[j])*cos(nuVec[j])*dnudphi0Vec[j];

		double dzOutdphi0 = a*cos(nuVec[j])*sinh(muVec[j])*dmudphi0Vec[j]
							- a*cosh(muVec[j])*sin(nuVec[j])*dnudphi0Vec[j];

		// except for the axis compute the values of the surface normals
		for(int i = 0; i < Nr; i++)
		{

			vur[i][j][0] = rVec[j]*cos(phi[i][j]);
			vur[i][j][1] = rVec[j]*sin(phi[i][j]);
			vur[i][j][2] = -zMean + zVec[j];

			vphi0[i][j][0] = ur[i][j]*drOutdphi0*cos(phi[i][j]) - r[i][j]*sin(phi[i][j])*dphidphi0Vec[j];
			vphi0[i][j][1] = ur[i][j]*drOutdphi0*sin(phi[i][j]) + r[i][j]*cos(phi[i][j])*dphidphi0Vec[j];
			vphi0[i][j][2] = ur[i][j]*dzOutdphi0;

			// first compute the vector in the normal direction which is not yet normalized
			crossProduct3( vur[i][j], vphi0[i][j], vCross[i][j] );

		}

	}

}



// requires that lid positions were already computed
void LidCDMLV::computeDisplacement()
{

	for(int i = 0; i < Nr; i++)
	{
		for(int j = 0; j < Nphi; j++)
		{

			uLid[i][j][0] = x[i][j] - x0[i][j];
			uLid[i][j][1] = y[i][j] - y0[i][j];
			uLid[i][j][2] = z[i][j] - z0[i][j];

		}
	}

}


void LidCDMLV::computeSecondSurfacePosition()
{

	// values at the endocardium

	for(int j = 0; j < Nphi; j++)
	{
		rVec[j] = a*sinh(muVec[j])*sin(nuVec[j]);
		zVec[j] = a*cosh(muVec[j])*cos(nuVec[j]);
	}


	// compute the mean value of z
	zMean = 0;
	for(int j = 0; j < Nphi; j++)
	{
		double cphi = computePeriodicIntegralCoef ( j, Nphi );
		zMean = zMean + cphi*zVec[j]*dphidphi0Vec[j]*dphi0;
	}
	zMean = zMean/(2*PI);


	for(int i = 0; i < Nr; i++)
	{

		for(int j = 0; j < Nphi; j++)
		{

			// r is computed to be the same as the value of r in the first surface
			r[i][j] = a*sinh(muSurf1[i][j])*sin(nuSurf1[i][j]);

			// ur coordinate from 0 to 1 must be calculated from the r grid of the first surface
			ur[i][j] = r[i][j]/rVec[j];

			// new z surface
			zLid[i][j] = (1-ur[i][j])*zMean + ur[i][j]*zVec[j];


			x[i][j] = r[i][j] * cos(phi[i][j]);
			y[i][j] = r[i][j] * sin(phi[i][j]);
			z[i][j] = zLid[i][j];

		}
	}

}




void LidCDMLV::computeFirstSurfacePosition( vector<double> & q, ParamCDMLV & prm )
{

	// compute the positions at the boundary
	for(int k = 0; k < Nphi; k++ )
	{

		double snu0 = sin(nu0Vec[k]);
		double cnu0 = cos(nu0Vec[k]);

		double cnu0sq = pow(cnu0,2);
		double snu0sq = pow(snu0,2);

		double cphi0 = cos(phi0Vec[k]);
		double sphi0 = sin(phi0Vec[k]);

		double fc0 = fc0Vec[k];
		double dfc0dnu0 = dfc0dnu0Vec[k];
		double dfc0dphi0 = dfc0dphi0Vec[k];

		double muin0 = mu0Vec[k]; // the vector mu0Vec stores values at the endocardium
		double nuUp0 = nuUp0Vec[k];

		double chmu0 = cosh(mu0Vec[k]);
		double shmu0 = sinh(mu0Vec[k]);
		double chmu0sq = pow(chmu0,2);

		double dmuin0dnu0 = dmuin0dnu0Vec[k];
		double dmuin0dphi0 = dmuin0dphi0Vec[k];

		// outputs (not all necessary, but required to have definitions)
		double mu, nu, phi, shmu, chmu, shmusq, chmusq, snu, cnu, snusq, cnusq,
		fc, dfcdnu0, dfcdphi0, K, sinRatio,  dnudnu0, dnudphi0, dphidnu0, dphidphi0,
		d2nudnu02, d2nudphi0dnu0,  d2nudphi02, d2phidnu02, d2phidphi02, d2phidphi0dnu0,
		dmudmu0, dmudnu0,  dmudphi0;


		// use the LV deformation function library to compute the full deformation

		prm.fourierDef.computeDeformation( muin0, nu0Vec[k], phi0Vec[k], q.data(), nuUp0, fc0, acube, chmu0, shmu0, chmu0sq, cnu0sq, cphi0, sphi0,
			snu0, snu0sq, cnu0, dfc0dnu0, dfc0dphi0, dmuin0dnu0, dmuin0dphi0,
			mu, nu, phi, shmu, chmu, shmusq, chmusq, snu, cnu, snusq, cnusq,
			fc, dfcdnu0, dfcdphi0, K, sinRatio,  dnudnu0, dnudphi0, dphidnu0, dphidphi0,
			d2nudnu02, d2nudphi0dnu0,  d2nudphi02, d2phidnu02, d2phidphi02, d2phidphi0dnu0,
			dmudmu0, dmudnu0,  dmudphi0 );


		muVec[k] = mu;
		nuVec[k] = nu;
		phiVec[k] = phi;

		dnudphi0Vec[k] = dnudphi0;
		dphidphi0Vec[k] = dphidphi0;

		dmudphi0Vec[k] = dmudphi0;

	}


	// compute the values on the interior of the surface
	double nuAxis = nuUp0Min;
	for(int i = 0; i < NmuLid; i++)
	{
		for(int j = 0; j < Nphi; j++)
		{
			muSurf1[i][j] = muVec[j]/mu0Vec[j]*mu0Surf1[i][j];
			nuSurf1[i][j] = nuAxis + ( nuVec[j] - nuAxis )* mu0Surf1[i][j]/mu0Vec[j];
			phi[i][j] = phiVec[j];

			zSurf1[i][j] = a*cosh(muSurf1[i][j])*cos(nuSurf1[i][j]);
		}
	}

}



void LidCDMLV::computeFirstIrregularVolume(){

	double dmu, dphi, dV, integrand;


	Virreg1 = 0;

	for(int i = 0; i < NmuLid; i++)
	{

		for(int j = 0; j < Nphi; j++)
		{

			// Volume is computed partially in the deformed, partially in the reference frame
			// The spacing is equal but possibly (likely) depends on phi
			dmu = muSurf1[1][j] - muSurf1[0][j];
			dphi = dphidphi0Vec[j]*dphi0;

			dV = dS0[i][j]*dmu*dphi;
			integrand = cos(nuVec[j])*sinh(muSurf1[i][j])*(cos(2*nuVec[j]) - 3*cosh(2*muSurf1[i][j]) - 2.0)
						- cos(nuSurf1[i][j])*sinh(muSurf1[i][j])*(cos(2*nuSurf1[i][j]) - 3*cosh(2*muSurf1[i][j]) - 2.0);

			Virreg1 = Virreg1 + integrand*dV;
		}
	}

	Virreg1 = Virreg1 * pow(a,3)/6.0;

}


void LidCDMLV::computeSecondIrregularVolume()
{


	// compute the integral in the r direction first
	// this integral has uneven spacing so it is computed
	// using the trapezoid rule with uneven spacing
	double * integrand_r = new double [Nr];
	double * Ir = new double [Nphi];

	for(int j = 0; j < Nphi; j++)
	{

		// compute the integrand
		for(int i = 0; i < Nr; i++)
		{
			integrand_r[i] = ur[i][j]*(zLid[i][j] - zSurf1[i][j])*pow(rVec[j],2);

		}

		// compute the integral
		Ir[j] = 0;
		for(int i = 0; i < Nr-1; i++)
		{
			double dur = ur[i+1][j] - ur[i][j];
			Ir[j] = Ir[j] + dur/2.0 * ( integrand_r[i+1] + integrand_r[i] );
		}

	}


	// compute the integral in the phi direction
	Virreg2 = 0;
	for(int j = 0; j < Nphi; j++)
	{
		double cphi = computePeriodicIntegralCoef ( j, Nphi );
		Virreg2 = Virreg2 + cphi*Ir[j]*dphidphi0Vec[j]*dphi0;
	}



	// clean up
	delete [] integrand_r;
	delete [] Ir;

}





