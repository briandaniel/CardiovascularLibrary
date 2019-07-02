/*
 * gdmlv_node.cpp
 *
 *  Created on: Sep 11, 2018
 *      Author: brian
 *
 *copied from
 *
 * node.cpp
 *
 *  Created on: Sep 9, 2017
 *      Author: brian
 */


#include "gdmlv_node.hpp"


void NodeGdmLV::initializeTensors( ParamGdmLV * prm )
{
	udisp = new double [3];
	vcrossSide = new double [3];
	vcrossTop = new double [3];

	tractionForcingSide = new double [3];
	tractionForcingTop = new double [3];
	bodyForcing = new double [3];
	nSide = new double [3];
	nTop = new double [3];
	divSigma_cart = new double [3];

	F = new double *[3];
	Finv = new double *[3];
	FT = new double *[3];
	C = new double *[3];
	Cinv = new double *[3];

	Q = new double *[3];
	QT = new double *[3];

	E = new double *[3];
	Efib = new double *[3];
	estrain = new double *[3];
	efib = new double *[3];
	Se = new double *[3];
	Sa0 = new double *[3];

	sigma_cart = new double *[3];
	sigma_prol = new double *[3];

	for(int k = 0; k < 3; k++)
	{
		F[k] = new double[3];
		Finv[k] = new double[3];
		FT[k] = new double[3];
		C[k] = new double[3];
		Cinv[k] = new double[3];

		Q[k] = new double[3];
		QT[k] = new double[3];

		E[k] = new double[3];
		Efib[k] = new double[3];
		estrain[k] = new double[3];
		efib[k] = new double[3];
		Se[k] = new double [3];
		Sa0[k] = new double [3];

		sigma_cart[k] = new double [3];
		sigma_prol[k] = new double [3];
	}

	for(int i = 0; i < 3; i++)
	{
		vcrossTop[i] = 0.0;
		vcrossSide[i] = 0.0;
		udisp[i] = 0.0;

		tractionForcingSide[i] = 0.0;
		tractionForcingTop[i] = 0.0;
		bodyForcing[i] = 0.0;
		nSide[i] = 0.0;
		nTop[i] = 0.0;
		divSigma_cart[i] = 0.0;

		for(int j = 0; j<3; j++)
		{
			F[i][j] = 0;
			Finv[i][j] = 0;
			FT[i][j] = 0;
			C[i][j] = 0;
			Cinv[i][j] = 0;

			Q[i][j] = 0;
			QT[i][j] = 0;

			E[i][j] = 0;
			Efib[i][j] = 0;
			estrain[i][j] = 0;
			efib[i][j] = 0;
			Se[i][j] = 0;

			sigma_cart[i][j] = 0.0;
			sigma_prol[i][j] = 0.0;
		}
	}

	// dq dependent tensors
	Nq = prm->Nq;
	dVIdq = new double [Nq];
	dEdq = new double ** [Nq];
	dEfibdq = new double ** [Nq];
	dudq = new double * [Nq];

	Svq = new double ** [Nq];
	Saq = new double ** [Nq];

	for(int k = 0; k < Nq; k++)
	{
		dVIdq[k] = 0.0;
		dudq[k] = new double [3];
		dEdq[k] = new double * [3];
		dEfibdq[k] = new double * [3];
		Svq[k] = new double * [3];
		Saq[k] = new double * [3];

		for( int i = 0; i < 3; i++)
		{

			dudq[k][i] = 0.0;
			dEdq[k][i] = new double [3];
			dEfibdq[k][i] = new double [3];

			Svq[k][i] = new double [3];
			Saq[k][i] = new double [3];

			for (int j = 0; j < 3; j++)
			{
				dEdq[k][i][j] = 0.0;
				dEfibdq[k][i][j] = 0.0;
				Svq[k][i][j] = 0.0;
				Saq[k][i][j] = 0.0;

			}
		}
	}


}


void NodeGdmLV::destroyTensors( ){

	// delete vectors
	if ( udisp != NULL ) { delete [] udisp; }
	if ( vcrossSide != NULL ) { delete [] vcrossSide; }
	if ( vcrossTop != NULL ) { delete [] vcrossTop; }

	if ( bodyForcing != NULL ) { delete [] bodyForcing; }
	if ( tractionForcingTop != NULL ) { delete [] tractionForcingTop; }
	if ( tractionForcingSide != NULL ) { delete [] tractionForcingSide; }
	if ( nTop != NULL ) { delete [] nTop; }
	if ( nSide != NULL ) { delete [] nSide; }
	if ( divSigma_cart != NULL ) { delete [] divSigma_cart; }

	// cout <<  "Cleaning up node " << globalID << endl;
	free33Matrix( F );
	free33Matrix( Finv );
	free33Matrix( FT );
	free33Matrix( C );
	free33Matrix( Cinv );

	free33Matrix( Q );
	free33Matrix( QT );

	free33Matrix( E );
	free33Matrix( Efib );
	free33Matrix( estrain );
	free33Matrix( efib );

	free33Matrix( Se );
	free33Matrix( sigma_cart );
	free33Matrix( sigma_prol );


	// free larger arrays
	if ( dEdq != NULL){
	for(int k = 0; k < Nq; k++)
	{
		free33Matrix( dEdq[k] );
	}
	delete [] dEdq;
	}

	if ( dEfibdq != NULL){
	for(int k = 0; k < Nq; k++)
	{
		free33Matrix( dEfibdq[k] );
	}
	delete [] dEfibdq;
	}


	if ( Svq!= NULL){
	for(int k = 0; k < Nq; k++)
	{
		free33Matrix( Svq[k] );
	}
	delete [] Svq;
	}

	if ( Saq != NULL){
	for(int k = 0; k < Nq; k++)
	{
		free33Matrix( Saq[k] );
	}
	delete [] Saq;
	}

	if ( dudq != NULL){
	for(int k = 0; k < Nq; k++)
	{
		delete [] dudq[k];
	}
	delete [] dudq;
	}



}

// Set the NodeGdmLV location
void NodeGdmLV::setLocation(int globalID_in, int localID_in, ParamGdmLV * prm){


	globalID = globalID_in;
	localID = localID_in;

	// location index
	imu = (int) floor( ( (double) globalID )/( prm->Nnu*prm->Nphi) );
	jnu = (int) floor( ( (double) ( globalID - imu*prm->Nnu*prm->Nphi))/( prm->Nphi) );
	kphi = globalID - imu*prm->Nnu*prm->Nphi - jnu*prm->Nphi;

	// Compute phi0
    double dphi0 = 2 * PI / prm->Nphi;
	phi0 = dphi0*kphi;

	// Compute nu0
	nuUp0 = computeNuUp0( phi0, prm->nuUp0Spline );
	double dnu0 = (PI - nuUp0 ) / (prm->Nnu- 1);
	nu0 = nuUp0 + dnu0*jnu;

	// Compute muin0, muout0
	prm->lvRef.computeRefSurfaceValuesOnly( nu0, phi0, muin0, muout0 );

	// The umu0 coordinate ranges from [0,1] between muin0 and muout0
	double dumu0 = ( 1.0 )/(prm->Nmu - 1.0);
	umu0 = imu*dumu0;

	// compute mu0
	mu0 = computeMu0( muin0, muout0, umu0 );


	// cout << "nu0 = " << nu0 << "  phi0 = " << phi0 << "  muin0 = " << muin0 << "  muout0 = " << muout0 << "  mu0 = " << mu0 << endl;

	// check if on this is in the RV surface region
	rvSurfNode = checkRVSurf( prm->rvBoundaryCoef, umu0, nu0, phi0 );

	// cout << "nu0 = " << nu0 << "  phi0 = " << phi0 << "  muin0 = " << muin0 << "  muout0 = " << muout0 << "  mu0 = " << mu0 << "  umu0 = " << umu0 << "  rvSurfNode = " << rvSurfNode << endl;

	// this requires surface computations if it is on the RV surface
	surfNode = rvSurfNode;

	// Set whether or not this is on the endocardial surface
	endoSurfNode = 0;

	if ( imu == 0 ){
		surfNode = 1;
		endoSurfNode = 1;


	}



	// Set the local copies of parameter variables as needed
	a = prm->a;

	// End diastolic stress initially set to zero
	eff_ed = 0;


}


void NodeGdmLV::setPrescribedLocation( double mu0In, double nu0In, double phi0In, int globalID_in, int localID_in, ParamGdmLV * prm){

	globalID = globalID_in;
	localID = localID_in;

	// values
	mu0 = mu0In;
	nu0 = nu0In;
	phi0 = phi0In;

	// location index
	imu = -1;
	jnu = -1;
	kphi = -1;

	// Compute nuUp0
	nuUp0 = computeNuUp0( phi0, prm->nuUp0Spline );

	// Compute muin0, muout0
	prm->lvRef.computeRefSurfaceValuesOnly( nu0, phi0, muin0, muout0 );


	// The umu0 coordinate ranges from [0,1] between muin0 and muout0
	umu0 = ( muout0 - mu0In )/( muout0 - muin0 );

	// check if on this is in the RV surface region
	rvSurfNode = checkRVSurf( prm->rvBoundaryCoef, umu0, nu0, phi0 );

	// this requires surface computations if it is on the RV surface
	surfNode = rvSurfNode;

	// Set whether or not this is on the endocardial surface
	endoSurfNode = 0;

	// epicardial surface point
	if( fabs(umu0 - 1) < 1e-12 )
	{
		surfNode = 1;
	}

	// endocardial surface point
	if ( umu0 < 1e-12 ){
		surfNode = 1;
		endoSurfNode = 1;
	}

	// Set the local copies of parameter variables as needed
	a = prm->a;

	// End diastolic stress initially set to zero
	eff_ed = 0;


}


// this is for plot exports
void NodeGdmLV::setFixedLocation( int imuIn, int jnuIn, int kphiIn, int Nmu, int Nnu, int Nphi, ParamGdmLV * prm){

	imu = imuIn;
	jnu = jnuIn;
	kphi = kphiIn;

	// Compute phi0 (adjusted so that the plotting is closed)
    double dphi0 = 2 * PI / ( Nphi - 1);
	phi0 = dphi0*kphi;

	// Compute nu0
	nuUp0 = computeNuUp0( phi0, prm->nuUp0Spline );
	double dnu0 = (PI - nuUp0 ) / (Nnu - 1);
	nu0 = nuUp0 + dnu0*jnu;

	// Compute muin0, muout0
	prm->lvRef.computeRefSurfaceValuesOnly( nu0, phi0, muin0, muout0 );

	// The umu0 coordinate ranges from [0,1] between muin0 and muout0
	double dumu0 = ( 1.0 )/( Nmu - 1.0);
	umu0 = imu*dumu0;

	// surface node ?
	endoSurfNode = 0;
	surfNode = 0;

	if ( umu0 < 1e-12 ){
		surfNode = 1;
		endoSurfNode = 1;
	}
	// epicardial surface node ?
	if( fabs(umu0 - 1) < 1e-12 )
	{
		surfNode = 1;
	}


	// compute mu0
	mu0 = computeMu0( muin0, muout0, umu0 );

	// Set the local copies of parameter variables as needed
	a = prm->a;


}




void NodeGdmLV::initialComputations( ParamGdmLV * prm )
{
	// set to zero in case they are used without being needed
	vCrossTopNorm = 0.0;
	vCrossSideNorm = 0.0;

	chmue = cosh(prm->mue);
	acube = pow(a,3);

	prolate2xyz( mu0, nu0, phi0, a, x0, y0, z0);

	snu0 = sin(nu0);
	cnu0 = cos(nu0);

	snu0sq = pow(snu0,2);
	cnu0sq = pow(cnu0,2);

	shmu0 = sinh(mu0);
	chmu0 = cosh(mu0);

	shmu0sq = pow(shmu0,2);
	chmu0sq = pow(chmu0,2);

	gmu0 = a*sqrt( shmu0sq + snu0sq );
	gnu0 = gmu0;
	gphi0 = a*shmu0*snu0;

	gmult0 = gmu0*gnu0*gphi0;

	chmuin0 = cosh(muin0);
	chmuin0sq = pow(chmuin0,2);
	shmuin0 = sinh(muin0);

	cphi0 = cos(phi0);
	sphi0 = sin(phi0);

	// evaluate spline surface derivatives used in additional computations (derivative of the constant part is 0)
	prm->lvRef.computeRefSurfaceValuesAndDerivatives( nu0, phi0,
			 muin0,  muout0,  dmuin0dnu0,  dmuin0dphi0,  dmuout0dnu0, dmuout0dphi0);

	computeFiberRotation( prm );

	// Compute the integration weights
	computeIntegrationWeights( prm );

	// compute the fc0 function and its derivatives
	fc0func( acube, chmuin0, chmuin0sq, shmuin0, dmuin0dnu0, dmuin0dphi0,
			 cnu0, cnu0sq, snu0, cphi0, fc0, dfc0dnu0, dfc0dphi0 );

}


void NodeGdmLV::computeIntegrationWeights( ParamGdmLV * prm )
{
	double dmu0 = (muout0 - muin0) / (prm->Nmu - 1);

	double dnu0 = (PI - nuUp0 ) / (prm->Nnu- 1);

	double dphi0 = 2 * PI / prm->Nphi;


	double cimu0 = computeIntegralCoef( imu, prm->Nmu );
	double cjnu0 = computeIntegralCoef( jnu, prm->Nnu );
	double ckphi0 = computePeriodicIntegralCoef ( kphi, prm->Nphi );

	// The weight for a volume integral = (1.)(2.)(3.) with
	// 1. equal spacing deltas
	// 2. integration rule weights
	// 3. scale factors
	IV0weight = (dmu0*dnu0*dphi0)*(cimu0*cjnu0*ckphi0)*(gmu0*gnu0*gphi0);


	// The weight for a ( nu0 x phi0 ) area integral = (1.)(2.) with
	// 1. equal spacing deltas
	// 2. integration rule weights
	// note that there are no scale factors, because the scale is determined by the deformed configuration as well (computed later)
	IA0weight = (dnu0*dphi0)*(cjnu0*ckphi0);


	// Same thing for the top integral
	IA0weight_top = (dmu0*dphi0)*(cimu0*ckphi0);


}





void NodeGdmLV::computeSubvalues( )
{

	prolate2xyz( mu, nu, phi, a, x, y, z);

	gmu = a*sqrt( shmusq + snusq );
	gnu = gmu;
	gphi = a*shmu*snu;

}









void NodeGdmLV::printNode( ofstream &fileID ){

	std::string structName;
	std::string key;
	structName = "nodeStruct.node" + to_string(globalID) + '.';

	// 1. indices
	key = "i";
	printMatlabVariableSimple( fileID, structName + key, imu );
	key = "j";
	printMatlabVariableSimple( fileID, structName + key, jnu );
	key = "k";
	printMatlabVariableSimple( fileID, structName + key, kphi );

	// 2. reference prolate coordinates
	key = "mu0";
	printMatlabVariableSimple( fileID, structName + key, mu0 );
	key = "nu0";
	printMatlabVariableSimple( fileID, structName + key, nu0 );
	key = "phi0";
	printMatlabVariableSimple( fileID, structName + key, phi0 );

	// 3. reference cartesian coordinates
	key = "x0";
	printMatlabVariableSimple( fileID, structName + key, x0 );
	key = "y0";
	printMatlabVariableSimple( fileID, structName + key, y0 );
	key = "z0";
	printMatlabVariableSimple( fileID, structName + key, z0 );

	// 4. deformed prolate coordinates
	key = "mu";
	printMatlabVariableSimple( fileID, structName + key, mu );
	key = "nu";
	printMatlabVariableSimple( fileID, structName + key, nu );
	key = "phi";
	printMatlabVariableSimple( fileID, structName + key, phi );

	// 5. deformed cartesian coordinates
	key = "x";
	printMatlabVariableSimple( fileID, structName + key, x );
	key = "y";
	printMatlabVariableSimple( fileID, structName + key, y );
	key = "z";
	printMatlabVariableSimple( fileID, structName + key, z );

	// 6. surface normals
	key = "nSide";
	printMatlabArraySimple( fileID, structName + key, nSide, 3);
	key = "nTop";
	printMatlabArraySimple( fileID, structName + key, nTop, 3);

	// 7. body forcing
	key = "bodyForcing";
	printMatlabArraySimple( fileID, structName + key, bodyForcing, 3);

	// 8. surface traction forcing
	key = "tractionForcingTop";
	printMatlabArraySimple( fileID, structName + key, tractionForcingTop, 3);
	key = "tractionForcingSide";
	printMatlabArraySimple( fileID, structName + key, tractionForcingSide, 3);


	// 9. cartesian representation of the divergence of the Pk1 stress
	key = "divSigma_cart";
	printMatlabArraySimple( fileID, structName + key, divSigma_cart, 3);


	// 10. displacement vector
	key = "udisp";
	printMatlabArraySimple( fileID, structName + key, udisp, 3);


	// 11. deformation gradient tensor
	key = "F";
	printMatlab2DArraySimple( fileID, structName + key, F, 3, 3);


	// 12. Cauchy stress
	computeCauchyStress(); // compute cauchy stress
	key = "sigma_prol";
	printMatlab2DArraySimple( fileID, structName + key, sigma_prol, 3, 3);
	key = "sigma_cart";
	printMatlab2DArraySimple( fileID, structName + key, sigma_cart, 3, 3);


	// 13. fiber stretch
	key = "lambda";
	printMatlabVariableSimple( fileID, structName + key, lambda );

	// 14. green strain
	key = "Eprol";
	printMatlab2DArraySimple( fileID, structName + key, E, 3, 3);


}


