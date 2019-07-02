/*
 * gdmlv_aggregate_computations.cpp
 *
 *  Created on: Dec 7, 2018
 *      Author: brian
 */



#include "gdmlv_aggregate_computations.hpp"



void computeAggregatesSeries( vector<vector<double>> & qStore, GdmLV & lv, vector <double> & VlvVec,
		vector <double> & shortAxisRadiusVec, vector <double> & longAxisLengthVec, vector <double> & meanTwistAngleVec )
{

	for(int i = 0; i < qStore.size(); i++ )
	{
		computeAggregates(qStore[i],lv,VlvVec[i],shortAxisRadiusVec[i],longAxisLengthVec[i],meanTwistAngleVec[i]);
	}


}



void computeAggregates( vector<double> & q, GdmLV & lv, double & Vlv, double & shortAxisMeanRadius, double & longAxisLength, double & meanTwistAngle )
{

	// 0. Compute lv at the current grid
	lv.computeLV(q.data(), 0, 0, 0);

	// 1. Take volume from standard computation
	Vlv = lv.Vlv;

	// 2. Compute long axis length
	longAxisLength = computeLongAxisLength( q, lv );

	// 3. Compute short axis mean radius
	shortAxisMeanRadius = computeShortAxisMeanRadius( q, lv );

	// 4. Compute mean twist
	meanTwistAngle = computeMeanTwistAngle( lv );

}


// Assumes LV has already been evaluated at correct q
double computeLongAxisLength(  vector<double> & q, GdmLV & lv  )
{

	// Point at base is taken from the lid
	double xMeanLid = 0;
	double yMeanLid = 0;
	double zMeanLid = lv.lid.zMean; // the mean z value of the lid gives the endpoint for the long axis length

	// Local node for computations
	NodeGdmLV node;
	// Initialize the nodes and the static values
	node.initializeTensors( & lv.prm );

	// 2. Compute apex endocardial location
	double nu0 = PI;
	double phi0 = 0;
	double muin0, muout0;
	lv.prm.lvRef.computeRefSurfaceValuesOnly( nu0, phi0, muin0, muout0 );

	node.setPrescribedLocation( muin0, nu0, phi0, 0, 0, &lv.prm );
	node.initialComputations( &lv.prm );
	node.computeDeformation( q.data(), lv.fourierDef, &lv.prm );
	node.computeSubvalues();

	double xApexEndo = node.x;
	double yApexEndo = node.y;
	double zApexEndo = node.z;


	double longAxisLength = sqrt( pow( xMeanLid - xApexEndo, 2) + pow( yMeanLid - yApexEndo, 2) + pow( zMeanLid - zApexEndo, 2) );

	return longAxisLength;

}


/*
double computeShortAxisMeanRadius( GdmLV & lv )
{

	// compute the mean value of the radius using the lid vectors
	double rMean = 0;
	for(int j = 0; j < lv.lid.Nphi; j++)
	{
		double cphi = computePeriodicIntegralCoef ( j, lv.lid.Nphi );
		rMean = rMean + cphi*lv.lid.rVec[j]*lv.lid.dphidphi0Vec[j]*lv.lid.dphi0;
	}
	rMean = rMean/(2*PI);

	return rMean;
}
*/


double computeShortAxisMeanRadius( vector<double> & q, GdmLV & lv )
{

	// Local node for computations
	NodeGdmLV node;
	// Initialize the nodes and the static values
	node.initializeTensors( & lv.prm );

	double muApexEndo, muApexEpi;
	lv.prm.lvRef.computeRefSurfaceValuesOnly( PI, 0, muApexEndo, muApexEpi );
	double muApex = muApexEpi;

	double xApex, yApex, zApex;
	prolate2xyz( muApex, PI, 0, lv.prm.a, xApex, yApex, zApex);

	double muBaseEndo, muBaseEpi;
	lv.prm.lvRef.computeRefSurfaceValuesOnly( lv.prm.nuUp0Min, 0, muBaseEndo, muBaseEpi );
	double muBase = muBaseEpi;

	double xBase, yBase, zBase;
	prolate2xyz( muBase, lv.prm.nuUp0Min, 0, lv.prm.a, xBase, yBase, zBase);

	double zMid = (zBase + zApex)/2;

	double mu0,nu0,phi0;
	xyz2prolate( 0, 0, zMid, lv.prm.a, mu0, nu0, phi0 );
	phi0 = 0;

	double rMean = 0;
	for(int kphi = 0; kphi < lv.prm.Nphi; kphi++)
	{
		double x,y,z;

		double dphi0 = 2 * PI / lv.prm.Nphi;
		phi0 = dphi0*kphi;

		double muMidEndo, muMidEpi;
		lv.prm.lvRef.computeRefSurfaceValuesOnly( nu0, phi0, muMidEndo, muMidEpi );
		mu0 = muMidEndo;

		node.setPrescribedLocation( mu0, nu0, phi0, 0, 0, &lv.prm );
		node.initialComputations( &lv.prm );
		node.computeDeformation( q.data(), lv.fourierDef, &lv.prm );
		node.computeSubvalues();


		rMean = rMean + sqrt(pow(node.x,2) + pow(node.y,2));
	}

	rMean = rMean/lv.prm.Nphi;

	return rMean;

}


double computeMeanTwistAngle( GdmLV & lv  )
{


	// compute the mean value of phi using the lid vectors
	double delphiMean = 0;
	for(int j = 0; j < lv.lid.Nphi; j++)
	{
		double cphi = computePeriodicIntegralCoef ( j, lv.lid.Nphi );
		delphiMean = delphiMean + cphi*(lv.lid.phiVec[j] - lv.lid.phi0Vec[j])*lv.lid.dphidphi0Vec[j]*lv.lid.dphi0;

	}
	delphiMean = delphiMean/(2*PI);


	// the mean twist is the mean phi value because the phi change is zero at the apex
	return delphiMean;

}








