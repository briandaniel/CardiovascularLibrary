/*
 * HalfLumpSimpleActivation_output.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: brian
 */



#include "HalfLumpSimpleActivation.hpp"





void outputModelResultHalfSimple( int Nsteps, ParamHalfSimple & lumpPrm, double * tStore, double ** wStore, double ** auxStore, double * wEnd )
{
	// Get the number of processes
	int Nprocs,procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	if ( procID == ROOT_ID )
	{

		ofstream outFile( lumpPrm.outputFile.c_str() , ios::out );

		printMatlabArraySimple( outFile, "wEnd", wEnd, lumpPrm.Nvars);

		// output time array
		printMatlabArraySimple( outFile, "t", tStore, Nsteps);

		// output variables
	    // W = [ V_{lv}, V_{la}, q_{mv}, q_{aov}, \zeta_{mv}, \zeta_{aov},  P_{sa}, P_{sp} ]
		printMatlabArraySimple( outFile, "Vlv", wStore[0], Nsteps);
		printMatlabArraySimple( outFile, "Vla", wStore[1], Nsteps);
		printMatlabArraySimple( outFile, "qmv", wStore[2], Nsteps);
		printMatlabArraySimple( outFile, "qaov", wStore[3], Nsteps);
		printMatlabArraySimple( outFile, "zetamv", wStore[4], Nsteps);
		printMatlabArraySimple( outFile, "zetaaov", wStore[5], Nsteps);
		printMatlabArraySimple( outFile, "Psa", wStore[6], Nsteps);
		printMatlabArraySimple( outFile, "Psp", wStore[7], Nsteps);

		// output auxiliary variables
		// auxVars = [Plv, Pla, Psv, At_ventricle, At_atrium, Aeff_mv, Aeff_aov]
		printMatlabArraySimple( outFile, "Plv", auxStore[0], Nsteps);
		printMatlabArraySimple( outFile, "Pla", auxStore[1], Nsteps);
		printMatlabArraySimple( outFile, "Psv", auxStore[2], Nsteps);
		printMatlabArraySimple( outFile, "At_ventricle", auxStore[3], Nsteps);
		printMatlabArraySimple( outFile, "At_atrium", auxStore[4], Nsteps);
		printMatlabArraySimple( outFile, "Aeff_mv", auxStore[5], Nsteps);
		printMatlabArraySimple( outFile, "Aeff_aov", auxStore[6], Nsteps);


		outFile.close();

	}

}



// print params
void ParamHalfSimple::printParams(){

	// Get the number of processes
	int Nprocs,procID;
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	if ( procID == ROOT_ID )
	{


		ofstream outFile( paramOutputFile.c_str() , ios::out );


		// Lump models
		printMatlabVariableSimple( outFile, "Csp", Csp);
		printMatlabVariableSimple( outFile, "Csa", Csa);

		printMatlabVariableSimple( outFile, "Rsp", Rsp);
		printMatlabVariableSimple( outFile, "Rsa", Rsa);

		printMatlabVariableSimple( outFile, "Ppv_fixed", Ppv_fixed);
		printMatlabVariableSimple( outFile, "Psv_fixed", Psv_fixed);
		printMatlabVariableSimple( outFile, "Rpv", Rpv);


		// Varying elastance models
		printMatlabVariableSimple( outFile, "Emax_lv", Emax_lv);
		printMatlabVariableSimple( outFile, "Emin_lv", Emin_lv);
		printMatlabVariableSimple( outFile, "Ks_lv", Ks_lv);
		printMatlabVariableSimple( outFile, "V0_lv", V0_lv);

		printMatlabVariableSimple( outFile, "Emax_la", Emax_la);
		printMatlabVariableSimple( outFile, "Emin_la", Emin_la);
		printMatlabVariableSimple( outFile, "Ks_la", Ks_la);
		printMatlabVariableSimple( outFile, "V0_la", V0_la);


		// activation model
		printMatlabVariableSimple( outFile, "Ta_ventricle", Ta_ventricle);
		printMatlabVariableSimple( outFile, "Ts_ventricle", Ts_ventricle);
		printMatlabVariableSimple( outFile, "kt_ventricle", kt_ventricle);

		printMatlabVariableSimple( outFile, "Ta_atrium", Ta_atrium);
		printMatlabVariableSimple( outFile, "Ts_atrium", Ts_atrium);
		printMatlabVariableSimple( outFile, "kt_atrium", kt_atrium);


		// Valve models
		printMatlabVariableSimple( outFile, "rho", rho);

		printMatlabVariableSimple( outFile, "Aeff_max_mv", Aeff_max_mv);
		printMatlabVariableSimple( outFile, "Aeff_min_mv", Aeff_min_mv);
		printMatlabVariableSimple( outFile, "leff_mv", leff_mv);
		printMatlabVariableSimple( outFile, "Kvo_mv", Kvo_mv);
		printMatlabVariableSimple( outFile, "Kvc_mv", Kvc_mv);
		printMatlabVariableSimple( outFile, "deltaP_open_mv", deltaP_open_mv);
		printMatlabVariableSimple( outFile, "deltaP_close_mv", deltaP_close_mv);

		printMatlabVariableSimple( outFile, "Aeff_max_aov", Aeff_max_aov);
		printMatlabVariableSimple( outFile, "Aeff_min_aov", Aeff_min_aov);
		printMatlabVariableSimple( outFile, "leff_aov", leff_aov);
		printMatlabVariableSimple( outFile, "Kvo_aov", Kvo_aov);
		printMatlabVariableSimple( outFile, "Kvc_aov", Kvc_aov);
		printMatlabVariableSimple( outFile, "deltaP_open_aov", deltaP_open_aov);
		printMatlabVariableSimple( outFile, "deltaP_close_aov", deltaP_close_aov);


		// Other
		printMatlabVariableSimple( outFile, "Tc", Tc);
		printMatlabVariableSimple( outFile, "Ncyc", Ncyc);
		printMatlabVariableSimple( outFile, "Nvars", Nvars);
		printMatlabVariableSimple( outFile, "dt", dt);


		outFile.close();
	}
}

