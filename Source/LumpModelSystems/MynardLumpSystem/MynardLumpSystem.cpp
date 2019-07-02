/*
 * MynardLumpSystem.cpp
 *
 *  Created on: Sep 13, 2018
 *      Author: brian
 */

#include "MynardLumpSystem.hpp"




void computeMynard0DSystem( ParamMynard0D & lumpPrm )
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
	int Nsteps = ceil( lumpPrm.Ncyc*lumpPrm.Tc/lumpPrm.dt );

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

	// set initial values
	setDefaultInitialValuesMynard0D( w );


	// compute initial values of aux variables
	computeMynard0DSystemTD( t, w, lumpPrm, dwdt, auxVars );


	// store initial values
	storeValuesMynard0DSystem( i, t, w, auxVars, tStore, wStore, auxStore, lumpPrm.Nvars, lumpPrm.NauxVars );






	//--------------------MAIN SOLUTION LOOP--------------------//
	if(procID == ROOT_ID)
		cout << endl << endl << "****************************** STARTING MAIN SOLUTION LOOP ******************************" << endl << endl;
	while ( t < (lumpPrm.Ncyc*lumpPrm.Tc - lumpPrm.dt) )
	{

		// Indicate at the command line current computation progress
		if( mod((double)i,100.0) == 0 )
		{
			if(procID == ROOT_ID)
			cout << "Computing steps " << i+1 << "-" << i+100 << " / " << Nsteps <<
					", simulation t = " << t << " sec / " << lumpPrm.Ncyc*lumpPrm.Tc << " sec. " << endl;
		}

		// update by euler step
		computeMynardTrapezoidIntegrationStep(t, w, lumpPrm, auxVars );

		t = t+lumpPrm.dt;
		i = i+1;


		// store current values
		storeValuesMynard0DSystem( i, t, w, auxVars, tStore, wStore, auxStore, lumpPrm.Nvars, lumpPrm.NauxVars );


	}
	//------------------END MAIN SOLUTION LOOP------------------//
	if(procID == ROOT_ID)
		cout << endl << "****************************** MAIN SOLUTION LOOP COMPLETED ******************************" << endl << endl << endl;
	print1DArrayLine( w, lumpPrm.Nvars, 8, "w");

	// Output parameters used in simulation
	lumpPrm.printParams();
	// Output result
	outputModelResultMynard0D( Nsteps, lumpPrm, tStore, wStore, auxStore, w );


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






void computeMynardTrapezoidIntegrationStep( double tn, double * w, ParamMynard0D & lumpPrm, double * auxVars )
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
	computeMynard0DSystemTD( tn, w, lumpPrm, Fn, auxVars );
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
		computeMynardGFunc( tnp1, lumpPrm.dt, w, X, Fn, lumpPrm, auxVars, GX );

		// compute the jacobian
		computeMynardApproximateJacobian( GX, X, lumpPrm.dX, J,
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



void computeMynardApproximateJacobian( double * GX, double * X, double dX, double ** J,
		double tnp1, double dt, double * wn, double * Fn, ParamMynard0D & lumpPrm, double * auxVars )
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
        computeMynardGFunc( tnp1, dt, wn, XdXk, Fn, lumpPrm, auxVars, GdXk );

    	for(int i = 0; i < lumpPrm.Nvars; i++)
    	{
    		J[i][k] = ( GdXk[i] - GX[i] )/dX;  // Compute approximate derivative
    	}

	}

    delete GdXk;
    delete [] XdXk;

}




void computeMynardGFunc( double tnp1, double dt, double * wn, double * wnp1_approx, double * Fn, ParamMynard0D & lumpPrm, double * auxVars, double * G )
{

	double * Fnp1_approx = new double [lumpPrm.Nvars];
	computeMynard0DSystemTD( tnp1, wnp1_approx, lumpPrm, Fnp1_approx, auxVars );

    // G is simply the finite difference approximation
	for(int k = 0; k < lumpPrm.Nvars; k++)
	{
		G[k] = wnp1_approx[k] - wn[k] - dt/2.0*( Fnp1_approx[k]  + Fn[k] );
	}

	// clean up
	delete [] Fnp1_approx;

}


void computeMynard0DSystemTD( double t, double * w, ParamMynard0D & lumpPrm, double * dwdt, double * auxVars )
{

	double Vlv, Vrv, Vla, Vra, qmv, qaov, qtri, qpul, zetamv, zetaaov, zetatri, zetapul, Ppv, Ppp, Ppa, Psp, Psa;
	double Plv, Prv, Pla, Pra, dVlvdt, dVrvdt, dVladt, dVradt;
	double qin_la, qin_ra;
	double At_ventricle, At_atrium;
	double Vu, VnotSV, Psv;
	double dqmv_dt, dZetamv_dt, dqaov_dt, dZetaaov_dt, dqtri_dt, dZetatri_dt, dqpul_dt, dZetapul_dt;

	//----------------- 0. pull variables -----------------//
    // W = [ Vlv, Vrv, Vla, Vra, qmv, qaov, qtri, qpul, zetamv, zetaaov, zetatri, zetapul, Ppv, Ppp, Ppa, Psp, Psa ]
	// Psv is not included as it is automatically determined by total blood volume conservation
	Vlv = w[0];
	Vrv = w[1];
	Vla = w[2];
	Vra = w[3];
	qmv = w[4];
	qaov = w[5];
	qtri = w[6];
	qpul = w[7];
	zetamv = w[8];
	zetaaov = w[9];
	zetatri = w[10];
	zetapul = w[11];
	Ppv = w[12];
	Ppp = w[13];
	Ppa = w[14];
	Psp = w[15];
	Psa = w[16];


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
	// 2. RV
	computeVaryingElastanceModelMynard( qtri, qpul, Vrv, At_ventricle,
			lumpPrm.Emax_rv, lumpPrm.Emin_rv, lumpPrm.Ks_rv, lumpPrm.V0_rv,
			Prv, dVrvdt );


	//----------------- 3. Compute reservoir pressure Psv to preserve total blood volume -----------------//
	// Zero pressure volume of the circulatory system
    Vu = lumpPrm.Vu_pv + lumpPrm.Vu_pp + lumpPrm.Vu_pa + lumpPrm.Vu_sv + lumpPrm.Vu_sp + lumpPrm.Vu_sa;
    // Volume of the rest of the circulation
    VnotSV = lumpPrm.Cpv*Ppv + lumpPrm.Cpp*Ppp + lumpPrm.Cpa*Ppa + lumpPrm.Csp*Psp + lumpPrm.Csa*Psa
    		+ Vlv + Vrv + Vla + Vra + Vu;
    // Pressure in system veins calculated to conserve total blood volume
    Psv = 1/lumpPrm.Csv*( lumpPrm.Vtotal - VnotSV );


	//----------------- 4. Compute varying elastance atrium models -----------------//

	// 1. LA
	computeVaryingElastanceModelAtrium0DCoupled( Ppv, lumpPrm.Rpv, qmv, Vla, At_atrium,
			lumpPrm.Emax_la, lumpPrm.Emin_la, lumpPrm.Ks_la, lumpPrm.V0_la,
			Pla, dVladt );

	// 2. RA
	computeVaryingElastanceModelAtrium0DCoupled( Psv, lumpPrm.Rsv, qtri, Vra, At_atrium,
			lumpPrm.Emax_ra, lumpPrm.Emin_ra, lumpPrm.Ks_ra, lumpPrm.V0_ra,
			Pra, dVradt );


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
	// 3. Tricuspid valve
	computeMynardValve( qtri, zetatri, Pra, Prv,
			lumpPrm.rho, lumpPrm.Aeff_max_tri, lumpPrm.Aeff_min_tri, lumpPrm.leff_tri,
			lumpPrm.Kvo_tri, lumpPrm.Kvc_tri, lumpPrm.deltaP_open_tri, lumpPrm.deltaP_close_tri,
			dqtri_dt, dZetatri_dt );
	// 4. Pulmonary valve
	computeMynardValve( qpul, zetapul, Prv, Ppa,
			lumpPrm.rho, lumpPrm.Aeff_max_pul, lumpPrm.Aeff_min_pul, lumpPrm.leff_pul,
			lumpPrm.Kvo_pul,	lumpPrm.Kvc_pul, lumpPrm.deltaP_open_pul, lumpPrm.deltaP_close_pul,
			dqpul_dt, dZetapul_dt );


	//----------------- 5. Set time derivatives -----------------//
    // W = [ Vlv, Vrv, Vla, Vra, qmv, qaov, qtri, qpul, zetamv, zetaaov, zetatri, zetapul, Ppv, Ppp, Ppa, Psp, Psa ]
	dwdt[0] = dVlvdt;
	dwdt[1] = dVrvdt;
	dwdt[2] = dVladt;
	dwdt[3] = dVradt;
	dwdt[4] = dqmv_dt;
	dwdt[5] = dqaov_dt;
	dwdt[6] = dqtri_dt;
	dwdt[7] = dqpul_dt;
	dwdt[8] = dZetamv_dt;
	dwdt[9] = dZetaaov_dt;
	dwdt[10] = dZetatri_dt;
	dwdt[11] = dZetapul_dt;


	//----------------- 6. Compute lumped parameter models -----------------//
    // W = [ Vlv, Vrv, Vla, Vra, qmv, qaov, qtri, qpul, zetamv, zetaaov, zetatri, zetapul, Ppv, Ppp, Ppa, Psp, Psa ]
	// 1. dPpvdt
	dwdt[12] = evaluateLumpModel( Ppv, Ppp, Pla, lumpPrm.Cpv, lumpPrm.Rpv, lumpPrm.Rpp);

	// 2. dPppdt
	dwdt[13] = evaluateLumpModel( Ppp, Ppa, Ppv, lumpPrm.Cpp, lumpPrm.Rpp, lumpPrm.Rpa);

	// 3. dPpadt
	dwdt[14] = 1/lumpPrm.Cpa*( qpul - (Ppa - Ppp)/lumpPrm.Rpa );

	// 4. dPspdt
	dwdt[15] = evaluateLumpModel( Psp, Psa, Psv, lumpPrm.Csp, lumpPrm.Rsp, lumpPrm.Rsa);

	// 5. dPsadt
	dwdt[16] = 1/lumpPrm.Csa*( qaov - (Psa - Psp)/lumpPrm.Rsa );


	//----------------- 7. Store aux variables -----------------//
	// auxVars = [Plv, Prv, Pla, Pra, Psv, At_ventricle, At_atrium]
	auxVars[0] = Plv;
	auxVars[1] = Prv;
	auxVars[2] = Pla;
	auxVars[3] = Pra;
	auxVars[4] = Psv;
	auxVars[5] = At_ventricle;
	auxVars[6] = At_atrium;


}







void storeValuesMynard0DSystem( int stepIdx, double t, double * w, double * auxVars,
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



