/*
 * GDMLVSystem_VRValves.cpp
 *
 *  Created on: Nov 20, 2018
 *      Author: brian
 */


#include "../GDMLVSystemVarResValves/GDMLVSystem_VRValves.hpp"




int computeSolution( GdmLV & lv, int procID, int Nprocs )
{

	int success = 0;

	// compute volume at zero pressure
	double Vlv0 = lv.computeZeroPressureVolume();

	// compute the maximum number of possible time steps
	int Nt = ( lv.prm.Ncyc*lv.prm.Tc ) / lv.prm.dtMin + 100;

	// storage arrays
	int Naux = 6;
	double ** qStore, ** wStore, *PlvStore, *PpaoStore, *VlvStore, *tStore, ** auxStore;


	// only use these on the root id
	if(procID == ROOT_ID)
	{
		qStore = new double * [ lv.prm.Nq ];
		for( int k = 0; k < lv.prm.Nq; k++ )
			qStore[k] = new double [ Nt ];

		wStore = new double * [ lv.prm.Nw ];
		for( int k = 0; k < lv.prm.Nw; k++ )
			wStore[k] = new double [ Nt ];

		auxStore = new double * [ Naux ];
		for( int k = 0; k < Naux; k++ )
			auxStore[k] = new double [ Nt ];

		PlvStore = new double [ Nt ];

		PpaoStore = new double [ Nt ];

		VlvStore = new double [ Nt ];

		tStore = new double [Nt];



		// set store variables to zero
		for(int k = 0; k < Nt; k++)
		{
			for( int j = 0; j < lv.prm.Nq; j++ )
				qStore[j][k] = 0.0;

			for( int j = 0; j < lv.prm.Nw; j++ )
				wStore[j][k] = 0.0;

			for( int j = 0; j < Naux; j++ )
				auxStore[j][k] = 0.0;


			PlvStore[k] = 0.0;
			PpaoStore[k] = 0.0;
			VlvStore[k] = 0.0;
			tStore[k] = 0.0;

		}
	}






	// Set up variables
	double Pla, Plv_guess, Ppao_guess, Plv, Ppao, Vlv;

	double * q = new double[ lv.prm.Nq ];
	double * dqdt = new double [lv.prm.Nq];
	double * w = new double [lv.prm.Nw];
	double * dwdt = new double [lv.prm.Nw];

	// set to initial values
	for(int k = 0; k < lv.prm.Nq; k++)
	{
		q[k] = lv.prm.q0[k];
		dqdt[k] = lv.prm.dqdt0[k];
	}

	for(int k = 0; k < lv.prm.Nw; k++)
	{
		w[k] = lv.prm.w0[k];
		dwdt[k] = lv.prm.dwdt0[k];
	}

	Plv = lv.prm.Plv0;
	Ppao = lv.prm.Ppao0;


	// --------------------------- main solution loop --------------------------- //
	double t = 0.0;
	int k = 0;
	double dt = 0;



	while( t <= lv.prm.Ncyc*lv.prm.Tc && Vlv == Vlv)
	{


		// compute time step
		dt = computeTimeStep( t, &lv.prm );

		// compute model
		// rk2Step(  q, w, t, dt, &lv, Plv, Ppao, Vlv, dqdt, dwdt );
		eulerStep(  q, w, t, dt, &lv, Plv, Ppao, Vlv, dqdt, dwdt );



		// output routine
		if(procID == ROOT_ID)
		{
			cout << "Computing step " << k << " / " << Nt << endl;
			cout << "  ";
			print1DArrayLine(q, lv.prm.Nq, 6, "q");
			cout << "  ";
			print1DArrayLine(w, lv.prm.Nw, 6, "w");
			cout << "  ";
			cout << "t = " << t << " sec,  " <<  "Plv = " << Plv <<  ", Vlv = " << Vlv  << ",  Ppao = " << Ppao <<  endl;
			cout << endl;


			for(int j = 0; j < lv.prm.Nq; j++)
			{
				qStore[j][k] = q[j];
			}
			for(int j = 0; j < lv.prm.Nw; j++)
			{
				wStore[j][k] = w[j];
			}

			PlvStore[k] = Plv;
			PpaoStore[k] = Ppao;
			VlvStore[k] = Vlv;
			tStore[k] = t;


			double Pla, Psa, Ppv, Ppa, Ppp, Psp, Pra, Vrv;
			double Psv, F_or, F_ir;
			double Rmv, Raov;
			pullW ( w, Pla, Psa, Ppv, Ppa, Ppp, Psp, Pra, Vrv );


			// this is needed to precompute the pressure in the right ventricle
			double Prv, Pmax_rv, F_or_dummy, F_ir_dummy, dVrv_dt_dummy;
			double At_rv = simpleSinActivation(t, lv.prm.Tc, lv.prm.Ta );
			computeVaryingElastanceModelUrsino( At_rv, Vrv, Pra, Ppa,
					lv.prm.Emax_rv, lv.prm.kr_rv, lv.prm.Rra, lv.prm.Vu_rv, lv.prm.ke_rv, lv.prm.P0_rv,
					Pmax_rv, Prv, F_ir_dummy, F_or_dummy, dVrv_dt_dummy);

			// compute Psv
			Psv = computePsv( Psa, Psp, Pra, Ppa, Ppp, Ppv, Pla, Vrv, Vlv, & lv.prm );

			// compute valve resistances
			valveResistances( Pla, Plv, Ppao, lv.prm.beta, lv.prm.Rmvc, lv.prm.Rmvo,
							lv.prm.Raovc, lv.prm.Raovo, Rmv, Raov );

			double At = simpleSinActivation(t, lv.prm.Tc, lv.prm.Ta );

			double AtLV = At;
			double AtRV = At;

			auxStore[0][k] = Rmv;
			auxStore[1][k] = Raov;
			auxStore[2][k] = AtLV;
			auxStore[3][k] = AtRV;
			auxStore[4][k] = Prv;
			auxStore[5][k] = Psv;

		}


		t = t+dt;
		k = k + 1;
	}


	int NtComputed = k-1;



	for(int i = 0; i < 3; i++)
	{
		lv.qFinal[i] = q[i];
	}


	// write the initial values
	if ( lv.prm.updateInitialValues == 1 )
	{
		printInitialValues ( procID, ROOT_ID, lv.prm.initialValuesFile, q, w,
						 dwdt, dqdt, Plv, Ppao, lv.prm.Nq, lv.prm.Nw );
	}


	// export variables on the root processor
	if(procID == ROOT_ID)
	{
		ofstream fileID( lv.prm.matOut.c_str() , ios::out );
		cout << "Writing output to matlab formatted file "  <<  lv.prm.matOut << " ... " << endl;

		printMatlabVariableSimple( fileID, "Nw", lv.prm.Nw);
		printMatlab2DArraySimple( fileID, "q", qStore, lv.prm.Nq, NtComputed);
		printMatlab2DArraySimple( fileID, "w", wStore, lv.prm.Nw, NtComputed );

		printMatlabArraySimple( fileID, "Plv", PlvStore, NtComputed );
		printMatlabArraySimple( fileID, "Ppao", PpaoStore, NtComputed );
		printMatlabArraySimple( fileID, "Vlv", VlvStore, NtComputed );
		printMatlabArraySimple( fileID, "t", tStore , NtComputed );
		printMatlab2DArraySimple( fileID, "auxVariables", auxStore, Naux, NtComputed );

		printMatlabVariableSimple( fileID, "Vlv0", Vlv0 );

		fileID.close();

		cout << "Finished writing." << endl;

		// print out the params
		lv.prm.printMatlabParamFile();


	}

	// clean up
	delete [] dqdt;
	delete [] dwdt;
	delete [] q;
	delete [] w;


	// only use these on the root id
	if(procID == ROOT_ID)
	{
		for( int k = 0; k < lv.prm.Nq; k++ )
			delete [] qStore[k];
		delete [] qStore;


		for( int k = 0; k < lv.prm.Nw; k++ )
			delete [] wStore[k];
		delete [] wStore;

		delete [] PlvStore;

		delete [] PpaoStore;

		delete [] tStore;

		delete [] VlvStore;

	}




	if ( t >= lv.prm.Ncyc*lv.prm.Tc )
		success = 1;

	return success;
}

void rk2Step( double * q, double * w, double t, double dt, GdmLV * lv,
					double &Plv, double &Ppao, double &Vlv, double *dqdt, double *dwdt )
{
	double * k1q = new double [lv->prm.Nq];
	double * k1w = new double [lv->prm.Nw];
	double * k2q = new double [lv->prm.Nq];
	double * k2w = new double [lv->prm.Nw];


	// set initial guess for k1 from previous time step
	for(int k = 0; k < lv->prm.Nq; k++)
		k1q[k] = dqdt[k];

	for(int k = 0; k < lv->prm.Nw; k++)
		k1w[k] = dwdt[k];


	evaluateModels( q, w, t, lv, Plv, Ppao, Vlv, k1q, k1w);

	// set intial guess for k2
	for(int k = 0; k < lv->prm.Nq; k++)
		k2q[k] = k1q[k];

	for(int k = 0; k < lv->prm.Nw; k++)
		k2w[k] = k2q[k];

	// update w/q per the rk2 method
	double * qdt2 = new double [lv->prm.Nq];
	double * wdt2 = new double [lv->prm.Nw];
	for(int k = 0; k < lv->prm.Nq; k++)
		qdt2[k] = q[k] + dt/2.0*k1q[k];

	for(int k = 0; k < lv->prm.Nw; k++)
		wdt2[k] = w[k] + dt/2.0*k1w[k];

	// compute k2
	evaluateModels( qdt2, wdt2, t + dt/2.0, lv, Plv, Ppao, Vlv, k2q, k2w);

	// compute result
	for(int k = 0; k < lv->prm.Nq; k++)
	{
		q[k] = q[k] + dt*k2q[k];
		dqdt[k] = k2q[k];
	}

	for(int k = 0; k < lv->prm.Nw; k++)
	{
		w[k] = w[k] + dt*k2w[k];
		dwdt[k] = k2w[k];

	}


	delete [] k1q;
	delete [] k1w;
	delete [] k2q;
	delete [] k2w;
	delete [] qdt2;
	delete [] wdt2;

}




void eulerStep( double * q, double * w, double t, double dt, GdmLV * lv,
					double &Plv, double &Ppao, double &Vlv, double *dqdt, double *dwdt )
{

	// evaluate
	evaluateModels( q, w, t, lv, Plv, Ppao, Vlv,  dqdt, dwdt);

	// update
	for(int k = 0; k < lv->prm.Nq; k++)
	{
		q[k] = q[k] + dt*dqdt[k];
	}

	for(int k = 0; k < lv->prm.Nw; k++)
	{
		w[k] = w[k] + dt*dwdt[k];
	}



}

void evaluateModels( double * q, double * w, double t, GdmLV * lv, double &Plv, double &Ppao, double &Vlv, double *dqdt, double *dwdt )
{

	double * dqdt_guess = new double [lv->prm.Nq];
	for(int k = 0; k < lv->prm.Nq; k++)
		dqdt_guess[k] = dqdt[k];

	double Plv_guess = Plv;
	double Ppao_guess = Ppao;

	double Pla, Psa, Ppv, Ppa, Ppp, Psp, Pra, Vrv;
	pullW ( w, Pla, Psa, Ppv, Ppa, Ppp, Psp, Pra, Vrv );

	// Activation function
	double At = simpleSinActivation(t, lv->prm.Tc, lv->prm.Ta  );
	double AtLV = At;

	// this is needed to precompute the pressure in the right ventricle
	double Prv, Pmax_rv, F_or_dummy, F_ir_dummy, dVrv_dt_dummy;

	double At_rv = At;
	computeVaryingElastanceModelUrsino( At_rv, Vrv, Pra, Ppa,
			lv->prm.Emax_rv, lv->prm.kr_rv, lv->prm.Rra, lv->prm.Vu_rv, lv->prm.ke_rv, lv->prm.P0_rv,
			Pmax_rv, Prv, F_ir_dummy, F_or_dummy, dVrv_dt_dummy);


	// update lv model

	lv->computeLV( q, t, AtLV, Prv );
	Vlv = lv->Vlv; // pull the volume



	if(lv->procID == ROOT_ID)
	{

		bool newtonSuccess;



		newtonSuccess = computeLVTimeDerivatives( q, t, Pla, Psa, lv->alpha, lv->eta,	lv->kappa, lv->dVlvdq,
									Plv_guess, Ppao_guess, dqdt_guess, &lv->prm, dqdt, Plv, Ppao);

		lumpTimeDerivatives( t, Plv, Ppao, Vlv, w, &lv->prm,  dwdt );


		/*
		cout << endl << endl << "--------------------------------------------------------" << endl;
		cout << "step t = " << t << endl;

		cout << "Ppao = " << Ppao << ", Plv = " << Plv << endl;

		print1DArrayLine( w, lv->prm.Nw, 12, "w" );
		print1DArrayLine( dwdt, lv->prm.Nw, 12, "dwdt" );

		print1DArrayLine( q, lv->prm.Nq, 12, "q" );
		print1DArrayLine( dqdt, lv->prm.Nq, 12, "dqdt" );

		cout << "--------------------------------------------------------" << endl << endl;
		*/

		//print1DArrayLine(lv->eta, lv->prm.Nq, 5, " eta" );
		//print1DArrayLine(lv->kappa, lv->prm.Nq, 5, " kappa" );

		double * F = new double [lv->prm.Nq];
		for(int i = 0; i < lv->prm.Nq; i++)
			F[i] = lv->kappa[i] - lv->eta[i]*Plv;
		//print1DArrayLine(F, lv->prm.Nq, 5, " F" );

		double normF = vectorNorm (F, lv->prm.Nq, 2);
		cout << "normF = " << normF << endl;

		delete [] F;


	}

	// broadcast the result to all processors
    MPI_Bcast(dqdt, lv->prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
    MPI_Bcast(q, lv->prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);

    MPI_Bcast(dwdt, lv->prm.Nw, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
    MPI_Bcast(w, lv->prm.Nw, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);

    MPI_Bcast(&Vlv, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
    MPI_Bcast(&Plv, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
    MPI_Bcast(&Ppao, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);


	delete [] dqdt_guess;
}







// This function assumes that the equilibrium integrals have already been evaluated
// inside the model
bool computeLVTimeDerivatives( double * q, double t, double Pla, double Psa, double ** alpha, double * eta,	double * kappa,
	double * dVlvdq, double Plv_guess, double Ppao_guess, double *dqdt_guess, ParamGdmLV * prm, double *dqdt, double &Plv, double &Ppao )
{

	double * Xprev = new double [prm->Nq+2];
	double * X = new double [prm->Nq+2];
	double * Xdiff = new double [prm->Nq+2];
	double * H = new double [prm->Nq+2];
	double ** J = new double * [prm->Nq+2];
	for(int i = 0; i < prm->Nq+2; i++)
		J[i] = new double [prm->Nq+2];

	// initial value of X
	for(int i = 0; i < prm->Nq; i++ )
		X[i] = dqdt_guess[i];


	X[prm->Nq] = Plv_guess;
	X[prm->Nq+1] = Ppao_guess;

	double Hnorm = prm->newtHConvg*2;
	double XdiffNorm = prm->newtXConvg*2;
	int k = 0;


	// compute the error
	lvNewtonJacobian( X, alpha, eta, dVlvdq,  Pla,Psa, prm, J);
	lvNewtonFunction( X, alpha, eta, kappa, dVlvdq, Pla, Psa, prm, H);
	// note: the LU solve screws up J during the solution
	luSolve( J, H, prm->Nq+2, Xdiff);

	XdiffNorm = vectorNorm ( Xdiff, prm->Nq+2, 2 );
	Hnorm = vectorNorm ( H, prm->Nq+2, 2 );

	// cout << endl << endl << "next newton iteration ..." << endl;
	while( (k < prm->newtMaxIter) && (Hnorm > prm->newtHConvg) && (XdiffNorm > prm->newtXConvg) )
	{

		// store X
		for(int i = 0; i < prm->Nq+2; i++ )
			Xprev[i] = X[i];

		// update X
		for(int i = 0; i < prm->Nq+2; i++ )
			X[i] = X[i] - Xdiff[i];

		if( Hnorm != Hnorm )
		{
			cout << "Newton iteration returned NaN. Solution may have failed. Using previous step in the iteration." << endl;

			for(int i = 0; i < prm->Nq+2; i++ )
				X[i] = Xprev[i];


			lvNewtonJacobian( X, alpha, eta, dVlvdq, Pla, Psa, prm, J);
			lvNewtonFunction( X, alpha, eta, kappa, dVlvdq, Pla, Psa, prm, H);

			// note: the LU solve screws up J during the solution
			luSolve( J, H, prm->Nq+2, Xdiff);

			print1DArrayLine( X, prm->Nq+2, 12, "X" );
			print1DArrayLine( H, prm->Nq+2, 12, "H" );


			exit(0);

			return 0;
		}


		// compute the jacobian/function values
		lvNewtonJacobian( X, alpha, eta, dVlvdq, Pla, Psa, prm, J);
		lvNewtonFunction( X, alpha, eta, kappa, dVlvdq, Pla, Psa, prm, H);

		// note: the LU solve screws up J during the solution
		luSolve( J, H, prm->Nq+2, Xdiff);

		// update the error norms
		XdiffNorm = vectorNorm ( Xdiff, prm->Nq+2, 2 );
		Hnorm = vectorNorm ( H, prm->Nq+2, 2 );


		k = k+1;
	}

	/*
	cout << endl << endl << "Newton solver ran " << k << " steps " << endl;
	cout << "output: " << endl;
	cout << "Pla = " << Pla << ", Psa = " << Psa << endl;

	print1DArrayLine( dVlvdq, prm->Nq, 12, "dVlvdq" );
	print1DArrayLine( X, prm->Nq+2, 12, "X" );
	print1DArrayLine( H, prm->Nq+2, 12, "H" );
	*/


	// save values
	for(int i = 0; i < prm->Nq; i++ )
		dqdt[i] = X[i];
	Plv = X[prm->Nq];
	Ppao = X[prm->Nq+1];


	// clean up
	delete [] X;
	delete [] Xprev;
	delete [] Xdiff;
	delete [] H;
	for(int i = 0; i < prm->Nq+2; i++)
		delete [] J[i];
	delete [] J;


	if (k >= prm->newtMaxIter)
		return 0;
	else
		return 1;

}


void lvNewtonJacobian( double * X, double ** alpha, double * eta,
		 	 	 	 	 double * dVlvdq, double Pla, double Psa, ParamGdmLV * prm, double ** J)
{
	double Plv, Ppao;
	double dHNq1_dPlv, dHNq1_dPpao, dHNq2_dPlv, dHNq2_dPpao;

	// pull values
	Plv = X[ prm->Nq ];
	Ppao = X[ prm->Nq + 1 ];

	// compute the derivatives of the nonlinear terms
	analyticNewtonDerivatives( Plv, Ppao, Pla, Psa, prm, dHNq1_dPlv, dHNq1_dPpao, dHNq2_dPlv, dHNq2_dPpao );

	for(int i = 0; i < prm->Nq+2; i++)
	{
		for(int j = 0; j < prm->Nq+2; j++)
		{
			J[i][j] = 0.0;
		}
	}

	for(int i = 0; i < prm->Nq; i++)
	{
		for(int j = 0; j < prm->Nq; j++)
		{
			J[i][j] = alpha[i][j];
		}
	}

	for( int i = 0; i < prm->Nq; i++)
	{
		J[i][prm->Nq] = -eta[i];
	}

	for( int j = 0; j < prm->Nq; j++)
	{
		J[prm->Nq][j] = dVlvdq[j];
	}


	J[prm->Nq][prm->Nq] = dHNq1_dPlv;
	J[prm->Nq][prm->Nq+1] = dHNq1_dPpao;
	J[prm->Nq+1][prm->Nq] = dHNq2_dPlv;
	J[prm->Nq+1][prm->Nq+1] = dHNq2_dPpao;


}



void analyticNewtonDerivatives( double Plv, double Ppao, double Pla, double Psa,
		ParamGdmLV * prm, double &dHNq1_dPlv, double &dHNq1_dPpao, double &dHNq2_dPlv, double &dHNq2_dPpao )
{
	double Rmv, Raov;
	double dRmv_dPlv, dRaov_dPlv, dRaov_dPpao;

	// compute valve resistances
	valveResistances( Pla, Plv, Ppao, prm->beta, prm->Rmvc, prm->Rmvo,
				prm->Raovc, prm->Raovo, Rmv, Raov );

	// analytic derivatives of the resistances
	valveResistanceDerivatives(  Pla, Plv, Ppao, prm->beta, prm->Rmvc, prm->Rmvo,
			prm->Raovc, prm->Raovo, dRmv_dPlv, dRaov_dPlv, dRaov_dPpao );

	// compute jacobian elements
	dHNq1_dPlv = ( Rmv + (Pla - Plv)*dRmv_dPlv )/( pow( Rmv,2) )
					+ ( Raov - (Plv - Ppao)*dRaov_dPlv )/( pow(Raov,2) );

	dHNq1_dPpao = - ( Raov + (Plv - Ppao)*dRaov_dPpao )/( pow(Raov,2) );

    dHNq2_dPlv = Ppao*dRaov_dPlv - prm->Rpao - Psa*dRaov_dPlv;

    dHNq2_dPpao = prm->Rpao + Raov + Ppao*dRaov_dPpao - Psa*dRaov_dPpao;

}




void lvNewtonFunction( double * X, double ** alpha, double * eta,
		double * kappa, double * dVlvdq, double Pla, double Psa, ParamGdmLV * prm, double * H ){

	double Plv, Ppao;
	double Rmv, Raov;

	// pull values
	Plv = X[ prm->Nq ];
	Ppao = X[ prm->Nq + 1 ];

	// compute valve resistances
	valveResistances( Pla, Plv, Ppao, prm->beta, prm->Rmvc, prm->Rmvo,
						prm->Raovc, prm->Raovo, Rmv, Raov );


	double * MX = new double [prm->Nq+2];
	double * rhs = new double [prm->Nq+2];
	double ** M = new double * [prm->Nq+2];
	for(int i = 0; i < prm->Nq+2; i++)
	{
		M[i] = new double [prm->Nq+2];
		for(int j = 0; j < prm->Nq+2; j++)
		{
			M[i][j] = 0.0;
		}
	}

	// Set up RHS
	for(int k = 0; k < prm->Nq; k++)
		rhs[k] = -kappa[k];
	rhs[prm->Nq] = Pla/Rmv;
	rhs[prm->Nq+1] = Psa*Raov;


	// set up the matrix
	for(int i = 0; i < prm->Nq; i++)
	{
		for(int j = 0; j < prm->Nq; j++)
		{
			M[i][j] = alpha[i][j];
		}
	}

	for( int i = 0; i < prm->Nq; i++)
	{
		M[i][prm->Nq] = -eta[i];
	}

	for( int j = 0; j < prm->Nq; j++)
	{
		M[prm->Nq][j] = dVlvdq[j];
	}

	M[prm->Nq][prm->Nq] = 1.0/Rmv + 1.0/Raov;
	M[prm->Nq][prm->Nq+1] = -1.0/Raov;
	M[prm->Nq+1][prm->Nq] = - prm->Rpao;
	M[prm->Nq+1][prm->Nq+1] = prm->Rpao + Raov;


	// compute MX
	matrixVectorMultiply( M, X, prm->Nq+2, prm->Nq+2, MX );

	for( int k = 0; k < prm->Nq + 2; k++)
	{
		H[k] = MX[k] - rhs[k];
	}



	// clean up
	delete [] MX;
	delete [] rhs;

	for(int k = 0; k < prm->Nq+2; k++)
	{
		delete [] M[k];
	}
	delete [] M;


}




double computeTimeStep( double t, ParamGdmLV * prm ){

	double dt;
	double c;

	double t_cycle = mod(t, prm->Tc);

	if( t_cycle < (prm->Tc-prm->Ta) )
	{
		dt = prm->dtMax;
	}else{

		dt = prm->dtMin;
	}

	return dt;

}



void pullW ( double * w, double & Pla, double & Psa, double & Ppv, double & Ppa, double & Ppp, double & Psp, double & Pra, double & Vrv )
{
	Pla = w[0];
	Psa = w[1];
	Ppv = w[2];
	Ppa = w[3];
	Ppp = w[4];
	Psp = w[5];
	Pra = w[6];
	Vrv = w[7];
}




void lumpTimeDerivatives( double t, double Plv, double Ppao, double Vlv, double* w, ParamGdmLV * prm,  double* dwdt )
{

		double Pla, Psa, Ppv, Ppa, Ppp, Psp, Pra, Vrv;
		double Prv, Psv, F_or, F_ir;
		double Rmv, Raov;
		double Pmax_rv, dVrv_dt;

		pullW ( w, Pla, Psa, Ppv, Ppa, Ppp, Psp, Pra, Vrv );

		// compute valve resistances
		valveResistances( Pla, Plv, Ppao, prm->beta, prm->Rmvc, prm->Rmvo,
				prm->Raovc, prm->Raovo, Rmv, Raov );

		// compute Psv
		Psv = computePsv ( Psa, Psp, Pra, Ppa, Ppp, Ppv, Pla, Vrv, Vlv, prm );
		Psv = 1.2;

		// compute right ventricle
		double At_rv = simpleSinActivation(t, prm->Tc, prm->Ta );
		computeVaryingElastanceModelUrsino( At_rv, Vrv, Pra, Ppa,
				prm->Emax_rv, prm->kr_rv, prm->Rra, prm->Vu_rv, prm->ke_rv, prm->P0_rv,
				Pmax_rv, Prv, F_ir, F_or, dVrv_dt);

		// compute time derivatives
		// d / dt[Pla]
		dwdt[0] = 1.0/prm->Cla*( (Ppv - Pla) / prm->Rpv - (Pla - Plv)/Rmv );

		// d / dt[Psa]
		dwdt[1] = 1.0 / prm->Csa*( (Ppao - Psa) / prm->Rpao  - (Psa - Psp)/prm->Rsa  );

		// d / dt[Ppv]
		dwdt[2] = 1.0 / prm->Cpv*((Ppp - Ppv) / prm->Rpp - (Ppv - Pla) / prm->Rpv);
		dwdt[2] = 0.0;
		w[2] = 1.2;

		// d / dt[Ppa]
		dwdt[3] = 1.0 / prm->Cpa*(F_or - (Ppa - Ppp)/prm->Rpa );

		// d / dt[Ppp]
		dwdt[4] = 1.0 / prm->Cpp*( (Ppa - Ppp)/prm->Rpa - (Ppp - Ppv) / prm->Rpp);

		// d / dt[Psp]
		dwdt[5] = 1.0 / prm->Csp*( (Psa - Psp)/prm->Rsa - (Psp - Psv) / prm->Rsp );

		// d / dt[Pra]
		dwdt[6] = 1.0 / prm->Cra * ((Psv - Pra) / prm->Rsv - F_ir);

		// d / dt[Vrv]
		dwdt[7] = F_ir - F_or;

}




// unstressed volume
double computePsv (  double Psa, double Psp, double Pra, double Ppa, double Ppp, double Ppv, double Pla, double Vrv, double Vlv, ParamGdmLV * prm )
{
	double Vu, Psv;

	Vu = prm->Vu_sa + prm->Vu_sp + prm->Vu_sv + prm->Vu_pa + prm->Vu_pp + prm->Vu_pv +  prm->Vu_ra + prm->Vu_la  ;

	Psv = 1.0 / prm->Csv*(prm->Vt0 - prm->Cla*Pla - Vlv - prm->Csa*Psa - prm->Csp*Psp
								   - prm->Cra*Pra - Vrv - prm->Cpa*Ppa - prm->Cpp*Ppp - prm->Cpv*Ppv - Vu);

	return Psv;
}




void printInitialValues ( int procID, int rootID, string initFileName, double * q, double * w,
							double * dwdt, double * dqdt, double Plv, double Ppao, int Nq, int Nw )
{


	// only print on rootID
	// and only if Plv is not NaN

	if((procID == rootID) && (Plv == Plv))
	{


		ofstream outFile(  initFileName.c_str() , ios::out );

		printMatlabVariableSimple( outFile, "Plv0", Plv );
		printMatlabVariableSimple( outFile, "Ppao0", Ppao );

		printMatlabArraySimple( outFile, "q0", q, Nq );
		printMatlabArraySimple( outFile, "w0", w, Nw );
		printMatlabArraySimple( outFile, "dqdt0", dqdt, Nq );
		printMatlabArraySimple( outFile, "dwdt0", dwdt, Nw );

		outFile.close();

		cout << "Exported last time step as new initial values." << endl;


	}
	else if ( Plv != Plv )
	{
		cout << "Values are NaN so the initial values were not reset." << endl;

	}


}



