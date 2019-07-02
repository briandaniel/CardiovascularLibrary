/*
 * SingleFiberMynardValveSystem.cpp
 *
 *  Created on: Dec 12, 2018
 *      Author: brian
 */


#include "SingleFiberMynardValveSystem.hpp"




void computeSingleFiberMynardValve( vector <double> & w, ParamSingleFiberMynardValve & lumpPrm )
{

	// Get the number of processes
	int Nprocs,procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	// preloop computations
	lumpPrm.paramPreloopComputations();

	// local variables
	vector <double> auxVars (lumpPrm.NauxVars,0);

	// number of steps
	int Nsteps = ceil( lumpPrm.Ncyc*lumpPrm.NtStepsCycle );

	vector<vector<double>> wStore ( lumpPrm.Nvars, vector<double>(Nsteps,0));
	vector<vector<double>> auxStore ( lumpPrm.NauxVars, vector<double>(Nsteps,0));
	vector<double> tStore (Nsteps,0);


	// compute solution
	computeSingleFiberMynardValveValues( w, lumpPrm, wStore, auxStore, tStore );


	// Output parameters used in simulation
	lumpPrm.printParams();
	// Output result
	outputModelResultSingleFiberMynardValve( Nsteps, wStore, auxStore, tStore, lumpPrm );


	if(getProcID() == ROOT_ID)
	{
		cout << "Solution completed with " << endl;
		cout << "w = "; print1DVector(w);
	}

}

void computeSingleFiberMynardValveValues( vector <double> & w, ParamSingleFiberMynardValve & lumpPrm,
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


	// store current values
	storeValuesSingleFiberMynardValve( w,  auxVars, t, i, lumpPrm, wStore, auxStore, tStore);


	// main solution loop
	while ( i < Nsteps)
	{

		// update w
		computeSingleFiberMynardValveTimeStep( t, dt, w, lumpPrm );


		t = t + dt;
		i = i+1;


		// store current values
		storeValuesSingleFiberMynardValve( w,  auxVars, t, i, lumpPrm, wStore, auxStore, tStore);

	}


}



void computeSingleFiberMynardValveTimeStep( double t, double dt, vector <double> & w,
		ParamSingleFiberMynardValve & lumpPrm )
{
	vector <double> wn(w.size(),0);
	vector <double> Fn(w.size()-2,0);
	for(int i = 0; i < w.size(); i++)
	{
		wn[i] = w[i];
	}
	double Gplvn, Gplan;
	computeSingleFiberMynardValveFFunc( wn, t, lumpPrm, Fn,	Gplvn, Gplan);



	// Object
	LumpSingleFiberMynardResidual Fobj;
	// 	void setPointers( double tIn, double dtIn, vector <double> & wn, vector <double> & Fn, ParamSingleFiberMynardValve & lumpPrm )
	Fobj.setPointers(t,dt,wn,Fn,lumpPrm);


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

	newt.findZero(w);


}




void computeSingleFiberMynardValveNewtonObj( vector<double> & wnp1, vector <double> & wn, vector <double > & Fn,
		double tn, double dt, ParamSingleFiberMynardValve lumpPrm, vector <double> & G )
{
	vector<double> Fnp1( lumpPrm.Nvars - 2, 0);
	double tnp1 = tn + dt;
	double Gplv, Gpla;

	computeSingleFiberMynardValveFFunc(  wnp1, tnp1, lumpPrm, Fnp1, Gplv, Gpla);

	for(int i = 0; i < lumpPrm.Nvars - 2; i++)
	{
		G[i] = wnp1[i] - wn[i] - dt/2.0*( Fnp1[i] + Fn[i] );
	}
	G[lumpPrm.Nvars-2] = Gplv;
	G[lumpPrm.Nvars-1] = Gpla;



}




void computeSingleFiberMynardValveFFunc( vector<double> & w, double t, ParamSingleFiberMynardValve lumpPrm, vector <double> & F,
		 double & Gplv, double & Gpla)
{


	double Vlv, Vla, qmv, qaov, zetamv, zetaaov, Psa, Psp;
	double Plv, Pla, dVlvdt, dVladt;
	double At_ventricle, At_atrium;
	double dqmv_dt, dZetamv_dt, dqaov_dt, dZetaaov_dt;
	double Ppv, Psv;
	double sigmaf_lv, sigmaf_la;
	double qinla;
	double PlvModel, PlaModel;

	//----------------- 0. pull variables -----------------//
	// W = [ V_{lv}, V_{la}, q_{mv}, q_{aov}, \zeta_{mv}, \zeta_{aov},  P_{sa}, P_{lv}, P_{la}  ]
	Vlv = w[0];
	Vla = w[1];
	qmv = w[2];
	qaov = w[3];
	zetamv = w[4];
	zetaaov = w[5];
	Psa = w[6];
	Plv = w[7];
	Pla = w[8];

	// Venous pressure set to fixed values
	Ppv = lumpPrm.Ppv_fixed;
	Psp = lumpPrm.Psp_fixed;

	//----------------- 1. Compute activation -----------------//
	// Ventricular activation
	At_ventricle = twoHillActivation( t, lumpPrm.m1_ventricle, lumpPrm.m2_ventricle, lumpPrm.tau1_ventricle,
			lumpPrm.tau2_ventricle, lumpPrm.Tc, lumpPrm.Ts_ventricle, lumpPrm.hillMaxVal_ventricle );
	// Atrial activation
	At_atrium = twoHillActivation( t, lumpPrm.m1_atrium, lumpPrm.m2_atrium, lumpPrm.tau1_atrium,
			lumpPrm.tau2_atrium, lumpPrm.Tc, lumpPrm.Ts_atrium, lumpPrm.hillMaxVal_atrium );


	//----------------- 2. Compute varying elastance ventricle models -----------------//
	// LV
	singleFiberModelGaussTension( qmv, qaov, Vlv, At_ventricle,	lumpPrm.V0_lv, lumpPrm.Vw_lv, lumpPrm.Ta0_lv, lumpPrm.Tp0_lv,
			lumpPrm.cp_lv, lumpPrm.cv_lv, lumpPrm.v0_lv, lumpPrm.ls0_lv, lumpPrm.lsmax_lv, lumpPrm.lsw_lv, PlvModel, dVlvdt, sigmaf_lv );

	// LA
	qinla = ( lumpPrm.Ppv_fixed - Pla )/ lumpPrm.Rpv;
	singleFiberModelGaussTension( qinla, qmv, Vla, At_atrium, lumpPrm.V0_la, lumpPrm.Vw_la, lumpPrm.Ta0_la, lumpPrm.Tp0_la,
			lumpPrm.cp_la, lumpPrm.cv_la, lumpPrm.v0_la, lumpPrm.ls0_la, lumpPrm.lsmax_la, lumpPrm.lsw_la, PlaModel, dVladt, sigmaf_la );

	//----------------- 5. Compute valve models -----------------//
	// 1. Mitral valve
	computeMynardValve( qmv, zetamv, Pla, Plv,
			lumpPrm.rho, lumpPrm.Aeff_max_mv, lumpPrm.Aeff_min_mv, lumpPrm.leff_mv,
			lumpPrm.Kvo_mv,	lumpPrm.Kvc_mv, lumpPrm.deltaP_open_mv, lumpPrm.deltaP_close_mv,
			dqmv_dt, dZetamv_dt );
	// 2. Aortic valve
	computeMynardValve( qaov, zetaaov, Plv, Psa,
			lumpPrm.rho, lumpPrm.Aeff_max_aov, lumpPrm.Aeff_min_aov, lumpPrm.leff_aov,
			lumpPrm.Kvo_aov,	lumpPrm.Kvc_aov, lumpPrm.deltaP_open_aov, lumpPrm.deltaP_close_aov,
			dqaov_dt, dZetaaov_dt );


	//----------------- 5. Set time derivatives -----------------//
	// W = [ V_{lv}, V_{la}, q_{mv}, q_{aov}, \zeta_{mv}, \zeta_{aov},  P_{sa}, P_{lv}, P_{la}  ]
	// F = dwdt
	F[0] = dVlvdt;
	F[1] = dVladt;
	F[2] = dqmv_dt;
	F[3] = dqaov_dt;
	F[4] = dZetamv_dt;
	F[5] = dZetaaov_dt;


	//----------------- 6. Compute lumped parameter models -----------------//
	// W = [ V_{lv}, V_{la}, q_{mv}, q_{aov}, \zeta_{mv}, \zeta_{aov},  P_{sa}, P_{lv}, P_{la}  ]
	// 1. dPsadt
	F[6] = 1/lumpPrm.Csa*( qaov - (Psa - Psp)/lumpPrm.Rsa );


	// These equations are not part of the ODE system
	// but are algebraic equations for the pressures that
	// have implicit definitions
	Gplv = Plv - PlvModel;
	Gpla = Pla - PlaModel;


}



void storeValuesSingleFiberMynardValve( vector<double> &w,  vector<double> &auxVars, double t, int stepIdx,
	ParamSingleFiberMynardValve & lumpPrm , vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector<double>& tStore)
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



















