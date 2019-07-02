/*
 * MynardHalfSystem.cpp
 *
 *  Created on: Oct 8, 2018
 *      Author: brian
 */




#include "MynardHalfSystem.hpp"




void computeMynardHalf0DSystem( double * w0, ParamMynardHalf0D & lumpPrm )
{

	// Get the number of processes
	int Nprocs,procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	// local variables
	double * w = new double [lumpPrm.Nvars];
	double * dwdt = new double [lumpPrm.Nvars];
	double * auxVars = new double [lumpPrm.NauxVars];

	// number of steps
	int Nsteps = ceil( lumpPrm.Ncyc*lumpPrm.NtStepsCycle );

	// Run preloop computations
	lumpPrm.paramPreloopComputations();

	// iteration variables
	double t = 0;
	int i = 0;

	// storage variables
	double * tStore = new double [Nsteps];
	double ** wStore = new double * [lumpPrm.Nvars];
	for(int k = 0; k < lumpPrm.Nvars; k++){	wStore[k] = new double [Nsteps]; }

	// aux storage
	double ** auxStore = new double * [lumpPrm.NauxVars];
	for(int k = 0; k < lumpPrm.NauxVars; k++){ auxStore[k] = new double [Nsteps]; }

	// set initial values to w0
	for(int k = 0; k < lumpPrm.Nvars; k++){ w[k] = w0[k]; }


	// compute initial values of aux variables
	computeMynardHalf0DSystemTD( t, w, lumpPrm, dwdt, auxVars );


	// store initial values
	storeValuesMynardHalf0DSystem( i, t, w, auxVars, tStore, wStore, auxStore, lumpPrm.Nvars, lumpPrm.NauxVars );






	//--------------------MAIN SOLUTION LOOP--------------------//
	if(procID == ROOT_ID)
		cout << endl << endl << "****************************** STARTING MAIN SOLUTION LOOP ******************************" << endl << endl;
	while ( i < lumpPrm.Ncyc*lumpPrm.NtStepsCycle )
	{

		// Indicate at the command line current computation progress
		if( mod((double)i,100.0) == 0 )
		{
			if(procID == ROOT_ID)
			cout << "Computing steps " << i+1 << "-" << i+100 << " / " << Nsteps <<
					", simulation t = " << t << " sec / " << lumpPrm.Ncyc*lumpPrm.Tc << " sec. " << endl;
		}

		// update by euler step
		computeMynardHalfTrapezoidIntegrationStep(t, w, lumpPrm, auxVars );

		t = t+lumpPrm.dt;
		i = i+1;


		// store current values
		storeValuesMynardHalf0DSystem( i, t, w, auxVars, tStore, wStore, auxStore, lumpPrm.Nvars, lumpPrm.NauxVars );


	}
	//------------------END MAIN SOLUTION LOOP------------------//
	if(procID == ROOT_ID)
	{
		cout << endl << "****************************** MAIN SOLUTION LOOP COMPLETED ******************************" << endl << endl << endl;
		print1DArrayLine( w, lumpPrm.Nvars, 8, "w");
	}
	// Output parameters used in simulation
	lumpPrm.printParams();
	// Output result
	outputModelResultMynardHalf0D( Nsteps, lumpPrm, tStore, wStore, auxStore, w );




	// clean up
	delete [] w;
	delete [] dwdt;
	delete [] auxVars;
	delete [] tStore;

	for(int k = 0; k < lumpPrm.Nvars; k++){	delete [] wStore[k]; }
	delete [] wStore;

	for(int k = 0; k < lumpPrm.NauxVars; k++){	delete [] auxStore[k]; }
	delete [] auxStore;

}




void computeMynardHalf0DSystemExport( double * w, double ** wStore, double ** auxStore, double * tStore, ParamMynardHalf0D & lumpPrm )
{

	// Get the number of processes
	int Nprocs,procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	// local variables
	double * dwdt = new double [lumpPrm.Nvars];
	double * auxVars = new double [lumpPrm.NauxVars];

	// number of steps
	int Nsteps = ceil( lumpPrm.Ncyc*lumpPrm.NtStepsCycle );

	// Run preloop computations
	lumpPrm.paramPreloopComputations();


	// iteration variables
	double t = 0;
	int i = 0;

	// compute initial values of aux variables
	computeMynardHalf0DSystemTD( t, w, lumpPrm, dwdt, auxVars );

	// store initial values
	storeValuesMynardHalf0DSystem( i, t, w, auxVars, tStore, wStore, auxStore, lumpPrm.Nvars, lumpPrm.NauxVars );


	//--------------------MAIN SOLUTION LOOP--------------------//

	while ( i < lumpPrm.Ncyc*lumpPrm.NtStepsCycle )
	{


		// update by euler step
		computeMynardHalfTrapezoidIntegrationStep(t, w, lumpPrm, auxVars );

		t = t+lumpPrm.dt;
		i = i+1;


		// store current values
		storeValuesMynardHalf0DSystem( i, t, w, auxVars, tStore, wStore, auxStore, lumpPrm.Nvars, lumpPrm.NauxVars );


	}
	//------------------END MAIN SOLUTION LOOP------------------//


	// clean up
	delete [] dwdt;
	delete [] auxVars;


}


void computeMynardHalfTrapezoidIntegrationStep( double tn, double * w, ParamMynardHalf0D & lumpPrm, double * auxVars )
{


	double tnp1 = tn + lumpPrm.dt;
	double * wnp1 = new double [lumpPrm.Nvars];
	double * GX = new double [lumpPrm.Nvars];
	double * X = new double [lumpPrm.Nvars];
	double * Xdiff = new double [lumpPrm.Nvars];
	double * Fn = new double [lumpPrm.Nvars];
	double ** J = new double * [lumpPrm.Nvars];
	for(int i = 0; i < lumpPrm.Nvars; i++)
	{
		J[i] = new double [lumpPrm.Nvars];
	}

	// compute starting guess using forward euler
	computeMynardHalf0DSystemTD( tn, w, lumpPrm, Fn, auxVars );
	for(int i = 0; i < lumpPrm.Nvars; i++)
	{
		X[i] = w[i] + lumpPrm.dt*Fn[i];
	}

	// newton iteration
	int k = 0;
	for(int i = 0; i < lumpPrm.Nvars; i++)
	{
	    Xdiff[i] = lumpPrm.xMinDiff*2;
		GX[i] = lumpPrm.fNewtMin*2;
	}

    while ( (vectorNorm(Xdiff,lumpPrm.Nvars,2.0) > lumpPrm.xMinDiff) &&  (vectorNorm(GX,lumpPrm.Nvars,2.0) > lumpPrm.fNewtMin) && (k <=  lumpPrm.maxNewtIter) )
	{

    	// compute the RHS
		computeMynardHalfGFunc( tnp1, lumpPrm.dt, w, X, Fn, lumpPrm, auxVars, GX );

		// compute the jacobian
		computeMynardHalfApproximateJacobian( GX, X, lumpPrm.dX, J,
				tnp1, lumpPrm.dt, w, Fn, lumpPrm, auxVars );

		// solve for the differences
		luSolve( J, GX, lumpPrm.Nvars, Xdiff );

		// update the newton iteration variable
    	for(int i = 0; i < lumpPrm.Nvars; i++)
    	{
    		X[i] = X[i] - Xdiff[i];
    	}

        k = k+1;

	}

	// update w
	for(int i = 0; i < lumpPrm.Nvars; i++)
	{
		w[i] = X[i];
	}

	// clean up
	delete [] wnp1;
	delete [] GX;
	delete [] X;
	delete [] Xdiff;
	delete [] Fn;
	for(int i = 0; i < lumpPrm.Nvars; i++)
	{
		delete [] J[i];
	}
	delete [] J;

}



void computeMynardHalfApproximateJacobian( double * GX, double * X, double dX, double ** J,
		double tnp1, double dt, double * wn, double * Fn, ParamMynardHalf0D & lumpPrm, double * auxVars )
{

	double * GdXk = new double [lumpPrm.Nvars];
	double * XdXk = new double [lumpPrm.Nvars];

    for ( int k = 0; k < lumpPrm.Nvars; k++ )
	{
    	for(int i = 0; i < lumpPrm.Nvars; i++)
    	{
    		XdXk[i] = X[i];
    	}
    	XdXk[k] = XdXk[k] + dX;

        // Xk is the approximation for wnp1
        // GdXk is the function to be solved for zero which is G
        computeMynardHalfGFunc( tnp1, dt, wn, XdXk, Fn, lumpPrm, auxVars, GdXk );

    	for(int i = 0; i < lumpPrm.Nvars; i++)
    	{
    		J[i][k] = ( GdXk[i] - GX[i] )/dX;  // Compute approximate derivative
    	}

	}

    delete GdXk;
    delete [] XdXk;

}




void computeMynardHalfGFunc( double tnp1, double dt, double * wn, double * wnp1_approx, double * Fn, ParamMynardHalf0D & lumpPrm, double * auxVars, double * G )
{

	double * Fnp1_approx = new double [lumpPrm.Nvars];
	computeMynardHalf0DSystemTD( tnp1, wnp1_approx, lumpPrm, Fnp1_approx, auxVars );

    // G is simply the finite difference approximation
	for(int k = 0; k < lumpPrm.Nvars; k++)
	{
		G[k] = wnp1_approx[k] - wn[k] - dt/2.0*( Fnp1_approx[k]  + Fn[k] );
	}

	// clean up
	delete [] Fnp1_approx;

}


void computeMynardHalf0DSystemTD( double t, double * w, ParamMynardHalf0D & lumpPrm, double * dwdt, double * auxVars )
{

	double Vlv, Vla, qmv, qaov, zetamv, zetaaov, Psa, Psp;
	double Plv, Pla, dVlvdt, dVladt;
	double At_ventricle, At_atrium;
	double dqmv_dt, dZetamv_dt, dqaov_dt, dZetaaov_dt;
	double Ppv, Psv;

	//----------------- 0. pull variables -----------------//
    // W = [ V_{lv}, V_{la}, q_{mv}, q_{aov}, \zeta_{mv}, \zeta_{aov},  P_{sa}, P_{sp} ]
	Vlv = w[0];
	Vla = w[1];
	qmv = w[2];
	qaov = w[3];
	zetamv = w[4];
	zetaaov = w[5];
	Psa = w[6];
	Psp = w[7];

	// Venous pressure set to fixed values
	Ppv = lumpPrm.Ppv_fixed;
	Psv = lumpPrm.Psv_fixed;


	//----------------- 1. Compute activation -----------------//
	// Ventricular activation
	At_ventricle = twoHillActivation( t, lumpPrm.m1_ventricle, lumpPrm.m2_ventricle, lumpPrm.tau1_ventricle,
			lumpPrm.tau2_ventricle, lumpPrm.Tc, lumpPrm.Ts_ventricle, lumpPrm.hillMaxVal_ventricle );
	// Atrial activation
	At_atrium = twoHillActivation( t, lumpPrm.m1_atrium, lumpPrm.m2_atrium, lumpPrm.tau1_atrium,
			lumpPrm.tau2_atrium, lumpPrm.Tc, lumpPrm.Ts_atrium, lumpPrm.hillMaxVal_atrium );


	//----------------- 2. Compute varying elastance ventricle models -----------------//
	// 1. LV
	computeVaryingElastanceModelMynard( qmv, qaov, Vlv, At_ventricle,
			lumpPrm.Emax_lv, lumpPrm.Emin_lv, lumpPrm.Ks_lv, lumpPrm.V0_lv,
			Plv, dVlvdt );


	//----------------- 3. Compute varying elastance atrium models -----------------//
	// 1. LA
	computeVaryingElastanceModelAtrium0DCoupled( Ppv, lumpPrm.Rpv, qmv, Vla, At_atrium,
			lumpPrm.Emax_la, lumpPrm.Emin_la, lumpPrm.Ks_la, lumpPrm.V0_la,
			Pla, dVladt );


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
    // W = [ V_{lv}, V_{la}, q_{mv}, q_{aov}, \zeta_{mv}, \zeta_{aov},  P_{sa}, P_{sp} ]
	dwdt[0] = dVlvdt;
	dwdt[1] = dVladt;
	dwdt[2] = dqmv_dt;
	dwdt[3] = dqaov_dt;
	dwdt[4] = dZetamv_dt;
	dwdt[5] = dZetaaov_dt;


	//----------------- 6. Compute lumped parameter models -----------------//
    // W = [ V_{lv}, V_{la}, q_{mv}, q_{aov}, \zeta_{mv}, \zeta_{aov},  P_{sa}, P_{sp} ]
	// 1. dPsadt
	dwdt[6] = 1/lumpPrm.Csa*( qaov - (Psa - Psp)/lumpPrm.Rsa );

	// 2. dPspdt
	dwdt[7] = evaluateLumpModel( Psp, Psa, Psv, lumpPrm.Csp, lumpPrm.Rsp, lumpPrm.Rsa);

	//----------------- 7. Store aux variables -----------------//
	// auxVars = [Plv, Pla, Psv, At_ventricle, At_atrium, Aeff_mv, Aeff_aov]
	auxVars[0] = Plv;
	auxVars[1] = Pla;
	auxVars[2] = Psv;
	auxVars[3] = At_ventricle;
	auxVars[4] = At_atrium;
	auxVars[5] = computeMynardEffectiveArea( zetamv, lumpPrm.Aeff_max_mv, lumpPrm.Aeff_min_mv );
	auxVars[6] = computeMynardEffectiveArea( zetaaov, lumpPrm.Aeff_max_aov, lumpPrm.Aeff_min_aov );


}







void storeValuesMynardHalf0DSystem( int stepIdx, double t, double * w, double * auxVars,
		double * tStore, double ** wStore, double ** auxStore, int Nvars, int NauxVars )
{

	// store time step
	tStore[stepIdx] = t;

	// store variables
	for(int k = 0; k < Nvars; k++)
	{
		wStore[k][stepIdx] = w[k];
	}

	// store aux variables
	for(int k = 0; k < NauxVars; k++)
	{
		auxStore[k][stepIdx] = auxVars[k];
	}

}
