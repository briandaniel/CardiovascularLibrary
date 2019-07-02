/*
 * CDMLV_StaticLoading.cpp
 *
 *  Created on: Jan 29, 2019
 *      Author: brian
 */


#include "CDMLV_StaticLoading.hpp"



// computes the parameters q that satisfy the variational equations under static loading conditions
void computeStaticSolution( CDMLVModel & lv, double Plv, double At, vector <double> & q, string prmPrefix, DataContainer & prmData )
{


	// Object
	StaticCDMLVResidual Fobj;
	Fobj.setPointers(lv, Plv, At ) ;

	// Set up variable
	int Nq = lv.prm.Nq;
	vector <double> X(Nq,0);
	for(int i = 0; i < Nq; i++)
	{
		X[i] = q[i];
	}

	// Algorithm
	NewtonsMethod newt;
	newt.setFobjPtr(Fobj);

	// void setParams( double xTolIn, double fTolIn, double dXJacobianIn, int maxIterIn, int verboseIn )
	double xTol = readPrmValue( "xTol_newt", prmPrefix, prmData);
	double fTol = readPrmValue( "fTol_newt", prmPrefix, prmData);
	double maxIter = readPrmValue( "maxIter_newt", prmPrefix, prmData);
	int verbose = (int) readPrmValue( "verbose_newt", prmPrefix, prmData);
	double dq = readPrmValue( "dq_newt", prmPrefix, prmData);

	newt.setParams( xTol, fTol, dq, maxIter, verbose);
	newt.findZero(X);


	for(int i = 0; i < Nq; i++)
	{
		q[i] = X[i];
	}


}



// computes the parameters q that satisfy the variational equations under static loading conditions
// computes incremental loading conditions from Plv = 0 up to Plv = Plv using steps deltaP
void computeStaticSolutionContinuation( double deltaP, CDMLVModel & lv, double Plv, double At, vector <double> & q, string prmPrefix, DataContainer & prmData )
{

	int Nq = lv.prm.Nq;
	double Pstep = 0;
	for(int i = 0; i < Nq; i++){ q[i] = 0.0; }

	int Nsteps = ceil(Plv/deltaP);
	double DeltaAt = At/Nsteps;
	int k = 0;

	while ( Pstep < Plv )
	{
		k = k+1;
		double AtStep = DeltaAt*k;

		if( Pstep + deltaP < Plv)
			Pstep = Pstep + deltaP;
		else
			Pstep = Plv;

		computeStaticSolution( lv, Pstep, AtStep, q, prmPrefix, prmData );

		if(q[0] != q[0])
			break;
	}


}


// computes the parameters q that satisfy the variational equations under static loading conditions
// computes incremental loading conditions from Plv = 0 up to Plv = Plv using steps deltaP
// Recomputes with small increments if sltn fails
void computeStaticSolutionContinuationRecomputes( double deltaP0,  CDMLVModel & lv, double Plv, double At, vector <double> & q, string prmPrefix, DataContainer & prmData )
{
	vector<double> qloc( q.size(), 0);
	for(int i = 0; i < q.size(); i++){ qloc[i] = q[i]; }

	double deltaP = deltaP0;
	computeStaticSolutionContinuation( deltaP, lv, Plv, At,  qloc, prmPrefix, prmData  );

	int n = 0;
	while(qloc[0] != qloc[0] && n < 10)
	{
		// reduce pressure step sizes in half
		deltaP = deltaP/2.0;

		// restore q to starting guess
		for(int i = 0; i < q.size(); i++){ qloc[i] = 0.0; }


		if(getProcID() == ROOT_ID)
			cout << "Recomputing static solution with smaller pressure steps to reach " << " Plv = " << Plv << endl;

		// try again
		computeStaticSolutionContinuation( deltaP, lv, Plv, At, qloc, prmPrefix, prmData );



		n++;
	}

	// return solution
	for(int i = 0; i < q.size(); i++){ q[i] = qloc[i]; }

}


