/*
 * CDMLV_HalfSystem.cpp
 *
 *  Created on: Jan 29, 2019
 *      Author: brian
 */

#include "CDMLV_HalfSystem.hpp"


void computeLVSystemHalf( CDMLVModel & lv, ParamLumpActive & lumpPrm, string initialValuesFileName, string outFileName )
{


	// Get the number of processes
	int Nprocs, procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	lumpPrm.paramPreloopComputations(lumpPrm.Tc);

	int Nw = lumpPrm.Nw;
	int Nq = lv.prm.Nq;
	int Naux = 11;


	// compute volume at zero pressure
	double Vlv0 = lv.computeZeroPressureVolume();

	// compute the maximum number of possible time steps
	int Nt = ( lumpPrm.Ncyc*lumpPrm.Tc ) / lumpPrm.dtMin + 100;


	vector<double> tStore (Nt,0);
	vector<vector<double>> varStore ( Nw+Nq+Naux, vector<double>(Nt,0));


	// storage arrays
	vector <double> X(Nq+3,0);
	vector <double> q(Nq,0);
	vector <double> dqdt(Nq,0);
	vector <double> w(Nw,0);
	vector <double> dwdt(Nw,0);

	// set to initial values
	readInitialValuesCDMLVHalf( initialValuesFileName, X, q, w  );


	if(procID == ROOT_ID )
	{
		cout << "Initial values: " << endl;
		cout << "   q = "; print1DVector(q);
		cout << "   w = "; print1DVector(w);
		cout << "   X = "; print1DVector(X);
		cout << endl;
	}

	// --------------------------- main solution loop --------------------------- //
	int NtComputed;
	computeLVSystemHalfVectorOutput( Naux, X, q, w, lv, lumpPrm, tStore, varStore, NtComputed );

	for(int i = 0; i < Nq; i++)
	{
		lv.qFinal[i] = q[i];
	}


	// export variables on the root processor
	exportDataValuesCDMLVHalf( outFileName, NtComputed, Nq, Nw, Naux, tStore, varStore );

	exportInitialValuesCDMLVHalf( initialValuesFileName, X, q, w );

	MPI_Barrier(MPI_COMM_WORLD);



}



void computeLVSystemHalfVectorOutput( int Naux, vector<double> & X, vector<double> & q, vector <double> & w,
		CDMLVModel & lv, ParamLumpActive & lumpPrm, vector<double> & tStore, vector<vector<double>> & varStore, int & NtComputed )
{
	int verbose = 1;


	int Nw = w.size();
	int Nq = q.size();
	vector <double> dqdt(Nq,0);
	vector <double> dwdt(Nw,0);
	vector <double> aux (Naux,0);

	if( getProcID() == ROOT_ID && verbose > 0 )
	{
		cout << "Initial values: " << endl;
		cout << "   q = "; print1DVector(q);
		cout << "   w = "; print1DVector(w);
		cout << "   X = "; print1DVector(X);
		cout << endl;
	}

	// --------------------------- main solution loop --------------------------- //
	double t = 0.0;
	int k = 0;
	double dt = 0;
	double Plv, Ppao, Pla, Vlv;

	evaluateModels_HalfSystem( t, q, w, lv, lumpPrm, X, Plv, Ppao, Pla, Vlv, dqdt, dwdt );
	saveDataValuesCDMLVHalf( 0, t, q, w,  aux, tStore, varStore, lumpPrm, Plv, Vlv, Pla, Ppao );


	while( t <= lumpPrm.Ncyc*lumpPrm.Tc && Vlv == Vlv)
	{


		// compute time step
		dt = computeTimeStepCDMLVHalf( t, lumpPrm );

		// compute model
		eulerStepCDMLVHalf( dt, t, q, w, lv, lumpPrm, X, Plv, Ppao, Pla, Vlv, dqdt, dwdt );

		// output routine
		if(getProcID() == ROOT_ID && mod(k,100) == 0 && verbose > 0 )
		{
			cout << "Computing time " << t << " / " << lumpPrm.Ncyc*lumpPrm.Tc << endl;
			cout << "   q = "; print1DVector(q);
			cout << "   w = "; print1DVector(w);
			cout << "dt = " << dt << " sec,  " <<  "Plv = " << Plv << ",  Pla = " << Pla <<
					", Vlv = " << Vlv  <<  ", Vla = " << w[0] << ",  Ppao = " << Ppao <<  endl;
			cout << endl;

		}

		t = t+dt;
		k = k + 1;

		saveDataValuesCDMLVHalf( k, t, q, w,  aux, tStore, varStore, lumpPrm, Plv, Vlv, Pla, Ppao );

	}

	NtComputed = k-1;


}

void eulerStepCDMLVHalf( double dt, double t, vector<double> & q, vector<double> & w, CDMLVModel & lv, ParamLumpActive & lumpPrm,
		vector <double> & X, double & Plv, double & Ppao, double & Pla, double & Vlv, vector <double> & dqdt,  vector <double> & dwdt )
{

	// evaluate
	evaluateModels_HalfSystem( t, q, w, lv, lumpPrm, X, Plv, Ppao, Pla, Vlv, dqdt, dwdt );

	// update
	for(int k = 0; k < lv.prm.Nq; k++)
	{
		q[k] = q[k] + dt*dqdt[k];
	}

	for(int k = 0; k < dwdt.size(); k++)
	{
		w[k] = w[k] + dt*dwdt[k];
	}



}


void evaluateModels_HalfSystem( double t, vector<double> & q, vector<double> & w, CDMLVModel & lv, ParamLumpActive & lumpPrm,
		vector <double> & X, double & Plv, double & Ppao, double & Pla, double & Vlv, vector <double> & dqdt,  vector <double> & dwdt )
{

	double Vla, Psa;
	// pull values
	// W = [ Vla, Psa ]
	Vla = w[0];
	Psa = w[1];


	// 1. Update LV model
	// Ventricular activation is instead a two-hill activation
	double At_ventricle = twoHillActivation( t, lumpPrm.m1_ventricle, lumpPrm.m2_ventricle, lumpPrm.tau1_ventricle + lumpPrm.activeDurationAdjustment,
			lumpPrm.tau2_ventricle - lumpPrm.activeDurationAdjustment, lumpPrm.Tc, lumpPrm.Ts_ventricle, lumpPrm.hillMaxVal_ventricle );
	double Prv = 0;


	lv.computeLV( q, t, At_ventricle );

	// no need to copy back to q as the values have not changed

	// 2. Compute time derivatives and pressures using a newton iteration
	computeCDMLVHalf_TimeDerivatives_Pressures( t, Vla, Psa, lumpPrm, lv, X, dqdt, Plv, Ppao, Pla );

	// 3. Update the lumped/single fiber model time derivatives
	lumpTimeDerivatives_CDMLVHalf( t, Plv, Ppao, Pla, lumpPrm, lv, w, dwdt );

	Vlv = lv.Vlv;

	// broadcast the result to all processors
    MPI_Bcast(X.data(), X.size(), MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);

    MPI_Bcast(dqdt.data(), dqdt.size(), MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
    MPI_Bcast(q.data(), q.size(), MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);

    MPI_Bcast(dwdt.data(), dwdt.size(), MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
    MPI_Bcast(w.data(), w.size(), MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);

    MPI_Bcast(&Vla, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
    MPI_Bcast(&Plv, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
    MPI_Bcast(&Ppao, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
    MPI_Bcast(&Vlv, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);

}


void lumpTimeDerivatives_CDMLVHalf( double t, double Plv, double Ppao, double Pla, ParamLumpActive & lumpPrm, CDMLVModel & lv,
		vector <double> & w, vector <double> & dwdt )
{
	double Vla, Psa;
	double At_atrium;
	double dVladt;
	double qmv, qaov;
	double Raov, Rmv;
	double qinla;
	double sigmaf_la;

	// pull values
	// W = [ Vla, Psa ]
	Vla = w[0];
	Psa = w[1];

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

	// LA
	qinla = ( lumpPrm.Ppv_fixed - Pla )/ lumpPrm.Rpv;
	singleFiberModelGaussTension( qinla, qmv, Vla, At_atrium,
			lumpPrm.V0_la, lumpPrm.Vw_la, lumpPrm.Ta0_la, lumpPrm.Tp0_la, lumpPrm.cp_la, lumpPrm.cv_la,
			lumpPrm.v0_la, lumpPrm.ls0_la, lumpPrm.lsmax_la, lumpPrm.lsw_la, Pla, dVladt, sigmaf_la );

	// W = [ Vla, Psa ]
	dwdt[0] = dVladt;
	dwdt[1] = 1/lumpPrm.Csa*(  (Ppao - Psa) / lumpPrm.Rpao  - (Psa - lumpPrm.Psp_fixed)/lumpPrm.Rsa );




}



void computeCDMLVHalf_TimeDerivatives_Pressures( double t, double Vla, double Psa, ParamLumpActive & lumpPrm, CDMLVModel & lv,
		vector <double> & X, vector <double> & dqdt, double & Plv, double & Ppao, double & Pla )
{
	int Nq = lv.prm.Nq;

	// Object
	CDMLVHalfNewtObj Fobj;
	Fobj.setPointers(t, Vla, Psa, lv, lumpPrm);

	vector<double> F(Nq+3,0);
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

	//Fobj.Feval(X,F);

	for(int i = 0; i < Nq; i++)
	{
		dqdt[i] = X[i];
	}
	Plv = X[Nq];
	Ppao = X[Nq+1];
	Pla = X[Nq+2];


}



void computeNewtonFunction( vector<double> & X, double t, double Vla, double Psa, CDMLVModel & lv, ParamLumpActive & lumpPrm, vector<double> & F)
{
	int Nq = lv.prm.Nq;
	double Plv, Ppao, Pla;
	double Rmv, Raov, qmv, qaov, qinla;
	double At_atrium, Pla_est, dVladt, sigmaf_la;
	vector<double> dqdt(Nq,0);

	// 0. Retrieve variables
	for(int i = 0; i < Nq; i++)
	{
		dqdt[i] = X[i];
	}
	Plv = X[Nq];
	Ppao = X[Nq+1];
	Pla = X[Nq+2];

	// 1. The first Nq equations are identical in form
	for(int i = 0; i < Nq; i++)
	{
		F[i] = lv.kappa[i] - lv.eta[i]*Plv;
		for(int j = 0; j < Nq; j++)
		{
			F[i] = F[i] + lv.alpha[i][j]*dqdt[j];
		}
	}

	// 2. The next equation is the volume conservation
	// valve resistances
	valveResistances( Pla, Plv, Ppao, lumpPrm.beta, lumpPrm.Rmvc, lumpPrm.Rmvo,
			lumpPrm.Raovc, lumpPrm.Raovo, Rmv, Raov );

	// estimated flows
	qmv = (Pla - Plv)/Rmv;
	qaov = (Plv - Ppao)/Raov;

	// compute the residual
	F[Nq] = -qmv + qaov;
	for(int j = 0; j < Nq; j++)
	{
		F[Nq] = F[Nq] + lv.dVlvdq[j]*dqdt[j];
	}

	// 3. The next equation is simply that for Ppao
	F[Nq+1] = Ppao*(lumpPrm.Rpao + Raov) - Plv*lumpPrm.Rpao - Psa*Raov;


	// 4. The final equation is built from the single fiber left atrium model
	// Atrial activation
	At_atrium = twoHillActivation( t, lumpPrm.m1_atrium, lumpPrm.m2_atrium, lumpPrm.tau1_atrium,
			lumpPrm.tau2_atrium, lumpPrm.Tc, lumpPrm.Ts_atrium, lumpPrm.hillMaxVal_atrium );

	qinla = ( lumpPrm.Ppv_fixed - Pla )/ lumpPrm.Rpv;
	singleFiberModelGaussTension( qinla, qmv, Vla, At_atrium,
			lumpPrm.V0_la, lumpPrm.Vw_la, lumpPrm.Ta0_la, lumpPrm.Tp0_la, lumpPrm.cp_la, lumpPrm.cv_la,
			lumpPrm.v0_la, lumpPrm.ls0_la, lumpPrm.lsmax_la, lumpPrm.lsw_la, Pla_est, dVladt, sigmaf_la );

	F[Nq+2] = Pla - Pla_est;

}



double computeTimeStepCDMLVHalf( double t, ParamLumpActive & lumpPrm ){

	double dt;
	double c;

	double t_cycle = mod(t, lumpPrm.Tc);

	if( t_cycle < (lumpPrm.Tc-lumpPrm.Ta) )
	{
		dt = lumpPrm.dtMax;
	}else{

		dt = lumpPrm.dtMin;
	}

	return dt;

}

