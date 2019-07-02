/*
 * SimpleLumpSystem.cpp
 *
 *  Created on: Sep 13, 2018
 *      Author: brian
 */

#include "SimpleLumpSystem.hpp"




void computeSimpleLumpSystem( ParamSimpleLump & lumpPrm )
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
	double * PlaStore = new double [Nsteps];
	double * PpvStore = new double [Nsteps];
	double * PppStore = new double [Nsteps];
	double * PpaStore = new double [Nsteps];
	double * PraStore = new double [Nsteps];
	double * PspStore = new double [Nsteps];
	double * PsaStore = new double [Nsteps];
	double * VlvStore = new double [Nsteps];
	double * VrvStore = new double [Nsteps];

	// aux storage
	double * PlvStore = new double [Nsteps];
	double * PrvStore = new double [Nsteps];
	double * Pmax_lvStore = new double [Nsteps];
	double * Pmax_rvStore = new double [Nsteps];
	double * PsvStore = new double [Nsteps];
	double * FilvStore = new double [Nsteps];
	double * FolvStore = new double [Nsteps];
	double * FirvStore = new double [Nsteps];
	double * ForvStore = new double [Nsteps];
	double * At_lvStore = new double [Nsteps];
	double * At_rvStore = new double [Nsteps];



	// set initial values
	setDefaultInitialValuesSimpleLump( w );

	// compute initial values of aux variables
	computeSimpleLumpSystemTD( t, w, lumpPrm, dwdt, auxVars );

	// store initial values
	storeValuesSimpleLump( i, t, w, auxVars, tStore, PlaStore, PpvStore,
			PppStore, PpaStore, PraStore, PspStore, PsaStore, VlvStore,
			VrvStore, PlvStore, PrvStore, Pmax_lvStore, Pmax_rvStore, PsvStore,
			FilvStore, FolvStore, FirvStore, ForvStore, At_lvStore, At_rvStore );




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
		computeSimpleLumpSystemTD( t, w, lumpPrm, dwdt, auxVars );
		for(int k = 0; k < lumpPrm.Nvars; k++)
		{
			w[k] = w[k] + lumpPrm.dt*dwdt[k];
		}

		t = t+lumpPrm.dt;
		i = i+1;


		// store current values
		storeValuesSimpleLump( i, t, w, auxVars, tStore, PlaStore, PpvStore,
				PppStore, PpaStore, PraStore, PspStore, PsaStore, VlvStore,
				VrvStore, PlvStore, PrvStore, Pmax_lvStore, Pmax_rvStore, PsvStore,
				FilvStore, FolvStore, FirvStore, ForvStore, At_lvStore, At_rvStore );


	}
	//------------------END MAIN SOLUTION LOOP------------------//
	if(procID == ROOT_ID)
		cout << endl << "****************************** MAIN SOLUTION LOOP COMPLETED ******************************" << endl << endl << endl;


	// Output parameters used in simulation
	lumpPrm.printParams();
	// Output result
    outputModelResultSimpleLump( Nsteps, lumpPrm, tStore, PlaStore, PpvStore,
			PppStore, PpaStore, PraStore, PspStore, PsaStore, VlvStore,
			VrvStore, PlvStore, PrvStore, Pmax_lvStore, Pmax_rvStore, PsvStore,
			FilvStore, FolvStore, FirvStore, ForvStore, At_lvStore, At_rvStore );


	// clean up
	delete [] w;
	delete [] dwdt;
	delete [] auxVars;
	delete [] tStore;
	delete [] PlaStore;
	delete [] PpvStore;
	delete [] PppStore;
	delete [] PpaStore;
	delete [] PraStore;
	delete [] PspStore;
	delete [] PsaStore;
	delete [] VlvStore;
	delete [] VrvStore;

	delete [] PlvStore;
	delete [] PrvStore;
	delete [] Pmax_lvStore;
	delete [] Pmax_rvStore;
	delete [] PsvStore;
	delete [] FilvStore;
	delete [] FolvStore;
	delete [] FirvStore;
	delete [] ForvStore;
	delete [] At_lvStore;
	delete [] At_rvStore;
}










void computeSimpleLumpSystemTD( double t, double * w, ParamSimpleLump & lumpPrm, double * dwdt, double * auxVars )
{

	double Pla, Ppv, Ppp, Ppa, Pra, Psv, Psp, Psa, Vlv, Vrv;
	double At, At_lv, At_rv;
	double Pmax_lv, Plv, Filv, Folv, dVlv_dt;
	double Pmax_rv, Prv, Firv, Forv, dVrv_dt;
	double Vu, VnotSV;

	//----------------- 0. pull variables -----------------//
    // W = [ Pla, Ppv, Ppp, Ppa, Pra, Psp, Psa, Vlv, Vrv ]
	// Psv is not included as it is automatically determined by total blood volume conservation
	Pla = w[0];
	Ppv = w[1];
	Ppp = w[2];
	Ppa = w[3];
	Pra = w[4];
	Psp = w[5];
	Psa = w[6];
	Vlv = w[7];
	Vrv = w[8];


	//----------------- 1. Compute activation -----------------//
	At = simpleSinActivation(t, lumpPrm.Tc, lumpPrm.Ta );
	At_lv = At;
	At_rv = At;


	//----------------- 2. Compute varying elastance models -----------------//
	// 1. LV
	computeVaryingElastanceModelUrsino( At_lv, Vlv, Pla, Psa,
			lumpPrm.Emax_lv, lumpPrm.kr_lv, lumpPrm.Rla, lumpPrm.Vu_lv, lumpPrm.ke_lv, lumpPrm.P0_lv,
			Pmax_lv, Plv, Filv, Folv, dVlv_dt);

	// 2. RV
	computeVaryingElastanceModelUrsino( At_rv, Vrv, Pra, Ppa,
			lumpPrm.Emax_rv, lumpPrm.kr_rv, lumpPrm.Rra, lumpPrm.Vu_rv, lumpPrm.ke_rv, lumpPrm.P0_rv,
			Pmax_rv, Prv, Firv, Forv, dVrv_dt);


	//----------------- 3. Compute reservoir pressure Psv to preserve total blood volume -----------------//
    Vu = lumpPrm.Vu_la + lumpPrm.Vu_pv + lumpPrm.Vu_pp + lumpPrm.Vu_pa + lumpPrm.Vu_ra
    		+ lumpPrm.Vu_sv + lumpPrm.Vu_sp + lumpPrm.Vu_sa + lumpPrm.Vu_lv + lumpPrm.Vu_rv;
    // Volume of the rest of the circulation
    VnotSV = lumpPrm.Cla*Pla + lumpPrm.Cpv*Ppv + lumpPrm.Cpp*Ppp + lumpPrm.Cpa*Ppa + lumpPrm.Cra*Pra
    		 + lumpPrm.Csp*Psp + lumpPrm.Csa*Psa + Vlv +  Vrv + Vu;
    // Pressure in system veins calculated to conserve total blood volume
    Psv = 1/lumpPrm.Csv*( lumpPrm.Vtotal - VnotSV );


	//----------------- 4. Set time derivatives -----------------//

	// 1. dPladt
	dwdt[0] = 1/lumpPrm.Cla*( (Ppv - Pla )/lumpPrm.Rpv - Filv );

	// 2. dPpvdt
	dwdt[1] = evaluateLumpModel( Ppv, Ppp, Pla, lumpPrm.Cpv, lumpPrm.Rpv, lumpPrm.Rpp);

	// 3. dPppdt
	dwdt[2] = evaluateLumpModel( Ppp, Ppa, Ppv, lumpPrm.Cpp, lumpPrm.Rpp, lumpPrm.Rpa);

	// 4. dPpadt
	dwdt[3] = 1/lumpPrm.Cpa*( Forv - (Ppa - Ppp)/lumpPrm.Rpa );

	// 5. dPradt
	dwdt[4] = 1/lumpPrm.Cra*( (Psv - Pra)/lumpPrm.Rsv - Firv );

	// 6. dPspdt
	dwdt[5] = evaluateLumpModel( Psp, Psa, Psv, lumpPrm.Csp, lumpPrm.Rsp, lumpPrm.Rsa);

	// 7. dPsadt
	dwdt[6] = 1/lumpPrm.Csa*( Folv - (Psa - Psp)/lumpPrm.Rsa );

	// 8-9. Previously computed dVlvdt and dVrvdt;
	dwdt[7] = dVlv_dt;
	dwdt[8] = dVrv_dt;



	//----------------- 5. Store aux variables -----------------//
	auxVars[0] = Plv;
	auxVars[1] = Prv;
	auxVars[2] = Pmax_lv;
	auxVars[3] = Pmax_rv;
	auxVars[4] = Psv;
	auxVars[5] = Filv;
	auxVars[6] = Folv;
	auxVars[7] = Firv;
	auxVars[8] = Forv;
	auxVars[9] = At_lv;
	auxVars[10] = At_rv;

}







void storeValuesSimpleLump( int stepIdx, double t, double * w, double * auxVars, double *	tStore, double * PlaStore, double * PpvStore,
		double * PppStore, double * PpaStore, double * PraStore, double * PspStore, double * PsaStore, double * VlvStore,
		double * VrvStore, double * PlvStore, double * PrvStore, double * Pmax_lvStore, double * Pmax_rvStore, double * PsvStore,
		double * FilvStore, double * FolvStore, double * FirvStore, double * ForvStore, double * At_lvStore, double * At_rvStore )
{

	// store variables
	tStore[stepIdx] = t;
	PlaStore[stepIdx] = w[0];
	PpvStore[stepIdx] = w[1];
	PppStore[stepIdx] = w[2];
	PpaStore[stepIdx] = w[3];
	PraStore[stepIdx] = w[4];
	PspStore[stepIdx] = w[5];
	PsaStore[stepIdx] = w[6];
	VlvStore[stepIdx] = w[7];
	VrvStore[stepIdx] = w[8];

	// store aux variables
	PlvStore[stepIdx] = auxVars[0];
	PrvStore[stepIdx] = auxVars[1];
	Pmax_lvStore[stepIdx] = auxVars[2];
	Pmax_rvStore[stepIdx] = auxVars[3];
	PsvStore[stepIdx] = auxVars[4];
	FilvStore[stepIdx] = auxVars[5];
	FolvStore[stepIdx] = auxVars[6];
	FirvStore[stepIdx] = auxVars[7];
	ForvStore[stepIdx] = auxVars[8];
	At_lvStore[stepIdx] = auxVars[9];
	At_rvStore[stepIdx] = auxVars[10];

}



