/*
 * SimpleLumpSystem.cpp
 *
 *  Created on: Sep 13, 2018
 *      Author: brian
 */

#include "SimpleLumpSystemHalf.hpp"


void computeSimpleLumpSystemHalf( vector <double> & w0, vector <double> & X0, ParamSimpleLumpHalf & lumpPrm )
{

	// Get the number of processes
	int Nprocs,procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	// preloop computations
	lumpPrm.paramPreloopComputations();

	// local variables
	vector<double> w (lumpPrm.Nvars,0);
	vector<double> dwdt (lumpPrm.Nvars,0);
	vector <double> auxVars (lumpPrm.NauxVars,0);
	vector <double> X(lumpPrm.Nnewt,0);

	// number of steps
	int Nsteps = ceil( lumpPrm.Ncyc*lumpPrm.NtStepsCycle );


	vector<vector<double>> wStore ( lumpPrm.Nvars, vector<double>(Nsteps,0));
	vector<vector<double>> auxStore ( lumpPrm.NauxVars, vector<double>(Nsteps,0));
	vector<double> tStore (Nsteps,0);


	// set initial values
	for(int k = 0; k < lumpPrm.Nvars; k++ )
	{
		w[k] = w0[k];
	}
	for( int k =0; k < lumpPrm.Nnewt; k++)
	{
		X[k] = X0[k];
	}


	// compute solution
	computeSimpleLumpSystemHalfValues( w, X, lumpPrm, wStore, auxStore, tStore );


	// Output parameters used in simulation
	lumpPrm.printParams();
	// Output result
	outputModelResultSimpleLumpHalf( Nsteps, wStore, auxStore, tStore, lumpPrm );


	// set initial values
	for(int k = 0; k < lumpPrm.Nvars; k++ )
	{
		w0[k] = w[k];
	}
	for( int k =0; k < lumpPrm.Nnewt; k++)
	{
		X0[k] = X[k];
	}

}


void computeSimpleLumpSystemHalfValues( vector <double> & w, vector <double> & X, ParamSimpleLumpHalf & lumpPrm,
		vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector <double> & tStore )
{

	// Get the number of processes
	int Nprocs,procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	// local variables
	vector<double> dwdt (lumpPrm.Nvars,0);
	vector <double> auxVars (lumpPrm.NauxVars,0);

	// number of steps
	int Nsteps = ceil( lumpPrm.Ncyc*lumpPrm.NtStepsCycle );
	double dt = lumpPrm.dt;

	// iteration variables
	double t = 0;
	int i = 0;

	// compute initial values of aux variables
	computeSimpleLumpSystemHalfTD( t, w, X, lumpPrm, dwdt, auxVars );
	storeValuesSimpleLumpHalf( w, auxVars, t, i, lumpPrm, wStore, auxStore, tStore);

	// main solution loop
	while ( i < Nsteps)
	{


		// Indicate at the command line current computation progress
		/*
		if( mod((double)i,100.0) == 0 )
		{
			if(procID == ROOT_ID)
			cout << "Computing steps " << i+1 << "-" << i+100 << " / " << Nsteps <<
					", simulation t = " << t << " sec / " << lumpPrm.Ncyc*lumpPrm.Tc << " sec. " << endl;
		}
		*/

		// update by euler step
		computeSimpleLumpSystemHalfTD( t, w, X, lumpPrm, dwdt, auxVars );
		for(int k = 0; k < lumpPrm.Nvars; k++)
		{
			w[k] = w[k] + dt*dwdt[k];
		}

		t = t + dt;
		i = i+1;


		// store current values
		storeValuesSimpleLumpHalf( w,  auxVars, t, i, lumpPrm, wStore, auxStore, tStore);

	}


}

void computeSimpleLumpSystemHalfTD( double t, vector <double> & w, vector<double> & X,
		ParamSimpleLumpHalf & lumpPrm, vector <double> & dwdt, vector<double> & auxVars )
{
	double Vla, Vlv, Psa;
	double At_ventricle;
	double At_atrium;
	double Plv, Pla, Ppao, dVlvdt, dVladt;
	double qmv, qaov;
	double Raov, Rmv;

	// pull values
    // W = [ Vla, Vlv, Psa ]
	Vla = w[0];
	Vlv = w[1];
	Psa = w[2];


	// Compute pressures via newton iteration
	computeSimpleLumpSystemPressures( t, w, X, lumpPrm, Plv, Pla, Ppao );


	// Ventricular activation
	At_ventricle = twoHillActivation( t, lumpPrm.m1_ventricle, lumpPrm.m2_ventricle, lumpPrm.tau1_ventricle,
			lumpPrm.tau2_ventricle, lumpPrm.Tc, lumpPrm.Ts_ventricle, lumpPrm.hillMaxVal_ventricle );
	// Atrial activation
	At_atrium = twoHillActivation( t, lumpPrm.m1_atrium, lumpPrm.m2_atrium, lumpPrm.tau1_atrium,
			lumpPrm.tau2_atrium, lumpPrm.Tc, lumpPrm.Ts_atrium, lumpPrm.hillMaxVal_atrium );

	// valve resistances
	valveResistances( Pla, Plv, Ppao, lumpPrm.beta, lumpPrm.Rmvc, lumpPrm.Rmvo,
			lumpPrm.Raovc, lumpPrm.Raovo, Rmv, Raov );

	// flows
	qmv = (Pla - Plv)/Rmv;
	qaov = (Plv - Ppao)/Raov;


	Ppao = (Plv*lumpPrm.Rpao + Psa*Raov)/(lumpPrm.Rpao + Raov);


	// LV
	// void computeVaryingElastanceModelMynard( double qin, double qout, double V, double At,
	//		double Emax, double Emin, double Ks, double V0, double & P, double & dVdt )
	computeVaryingElastanceModelMynard( qmv, qaov, Vlv, At_ventricle,
			lumpPrm.Emax_lv, lumpPrm.Emin_lv, lumpPrm.Ks_lv, lumpPrm.V0_lv, Plv, dVlvdt );

	// LA
	computeVaryingElastanceModelAtrium0DCoupled( lumpPrm.Ppv_fixed, lumpPrm.Rpv, qmv, Vla, At_atrium,
			lumpPrm.Emax_la, lumpPrm.Emin_la, lumpPrm.Ks_la, lumpPrm.V0_la,
			Pla, dVladt );


    // W = [ Vla, Vlv, Psa ]
	dwdt[0] = dVladt;
	dwdt[1] = dVlvdt;
	dwdt[2] = 1/lumpPrm.Csa*(  (Ppao - Psa) / lumpPrm.Rpao  - (Psa - lumpPrm.Psp_fixed)/lumpPrm.Rsa );


	//----------------- Store aux variables -----------------//
	// auxVars = [Plv, Pla, At_ventricle, At_atrium, qmv, qaov, Rmv, Raov ]
	auxVars[0] = Plv;
	auxVars[1] = Pla;
	auxVars[2] = At_ventricle;
	auxVars[3] = At_atrium;
	auxVars[4] = qmv;
	auxVars[5] = qaov;
	auxVars[6] = Rmv;
	auxVars[7] = Raov;



}



void computeSimpleLumpSystemPressures( double t, vector <double> & w, vector<double> & X,
		ParamSimpleLumpHalf & lumpPrm, double & Plv, double & Pla, double & Ppao )
{
	// Object
	LumpResidual Fobj;
	Fobj.setPointers(t,w,lumpPrm);

	vector<double> F(3,0);
	Fobj.Feval(X,F);

	// Algorithm
	NewtonsMethod newt;
	newt.setFobjPtr(Fobj);

	// void setParams( double xTolIn, double fTolIn, double dXJacobianIn, int maxIterIn, int verboseIn )
	string prmPrefix = "staticLV";
	double xTol = lumpPrm.xMinDiff;
	double fTol = lumpPrm.fNewtMin;
	double maxIter = lumpPrm.maxNewtIter;
	int verbose = 0;
	double dX = lumpPrm.dX;

	newt.setParams( xTol, fTol, dX, maxIter, verbose);

	newt.findZero(X);

	Fobj.Feval(X,F);

	Pla = X[0];
	Plv = X[1];
	Ppao = X[2];

}




void computeSimpleLumpSystemHalfNewtonObj( vector<double> &X, vector<double> &F,
		double t, vector <double> & w, ParamSimpleLumpHalf & lumpPrm )
{

	double Vla, Vlv, Psa;
	double At_ventricle;
	double At_atrium;
	double Plv, Pla, Ppao, dVlvdt, dVladt;
	double Plv_est, Pla_est, Ppao_est;
	double qmv, qaov;
	double Raov, Rmv;

	//----------------- 0. pull variables -----------------//
    // W = [ Vla, Vlv, Psa ]
	Vla = w[0];
	Vlv = w[1];
	Psa = w[2];

	// X = [Pla, Plv, Ppao]
	Pla = X[0];
	Plv = X[1];
	Ppao = X[2];

	// Ventricular activation
	At_ventricle = twoHillActivation( t, lumpPrm.m1_ventricle, lumpPrm.m2_ventricle, lumpPrm.tau1_ventricle,
			lumpPrm.tau2_ventricle, lumpPrm.Tc, lumpPrm.Ts_ventricle, lumpPrm.hillMaxVal_ventricle );
	// Atrial activation
	At_atrium = twoHillActivation( t, lumpPrm.m1_atrium, lumpPrm.m2_atrium, lumpPrm.tau1_atrium,
			lumpPrm.tau2_atrium, lumpPrm.Tc, lumpPrm.Ts_atrium, lumpPrm.hillMaxVal_atrium );

	// valve resistances
	valveResistances( Pla, Plv, Ppao, lumpPrm.beta, lumpPrm.Rmvc, lumpPrm.Rmvo,
			lumpPrm.Raovc, lumpPrm.Raovo, Rmv, Raov );

	// flows
	qmv = (Pla - Plv)/Rmv;
	qaov = (Plv - Ppao)/Raov;

	// What the model computes Ppao to be
	Ppao_est = (Plv*lumpPrm.Rpao + Psa*Raov)/(lumpPrm.Rpao + Raov);

	// LV
	// void computeVaryingElastanceModelMynard( double qin, double qout, double V, double At,
	//		double Emax, double Emin, double Ks, double V0, double & P, double & dVdt )
	computeVaryingElastanceModelMynard( qmv, qaov, Vlv, At_ventricle,
			lumpPrm.Emax_lv, lumpPrm.Emin_lv, lumpPrm.Ks_lv, lumpPrm.V0_lv, Plv_est, dVlvdt );

	// LA
	computeVaryingElastanceModelAtrium0DCoupled( lumpPrm.Ppv_fixed, lumpPrm.Rpv, qmv, Vla, At_atrium,
			lumpPrm.Emax_la, lumpPrm.Emin_la, lumpPrm.Ks_la, lumpPrm.V0_la,
			Pla_est, dVladt );

	// Estimate the function values: difference between current estimate and approximated
	F[0] = Pla - Pla_est;
	F[1] = Plv - Plv_est;
	F[2] = Ppao - Ppao_est;

}




void storeValuesSimpleLumpHalf( vector<double> &w,  vector<double> &auxVars, double t, int stepIdx,
		ParamSimpleLumpHalf & lumpPrm , vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector<double>& tStore)
{

		for(int i = 0; i < lumpPrm.Nvars; i++)
		{
			wStore[i][stepIdx] = w[i];
		}

		for(int i = 0; i < lumpPrm.NauxVars; i++)
		{
			auxStore[i][stepIdx] = auxVars[i];
		}
		tStore[stepIdx] = t;
}













