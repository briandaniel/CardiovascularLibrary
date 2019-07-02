/*
 * SingleFiberMynardValveSystem_output.cpp
 *
 *  Created on: Dec 12, 2018
 *      Author: brian
 */


#include "SingleFiberMynardValveSystem.hpp"

void outputModelResultSingleFiberMynardValve( int Nsteps, vector<vector<double>> & wStore,
		vector<vector<double>> & auxStore, vector<double>& tStore, ParamSingleFiberMynardValve & lumpPrm )
{

	if( getProcID() == 0)
	{
		int idx;

		ofstream outFile( lumpPrm.outputFile.c_str() , ios::out );

		// W = [ V_{lv}, V_{la}, q_{mv}, q_{aov}, \zeta_{mv}, \zeta_{aov},  P_{sa}, P_{lv}, P_{la}  ]
		idx = 0;
		printMatlabArraySimple( outFile, "Vlv", wStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "Vla", wStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "qmv", wStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "qaov", wStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "zetamv", wStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "zetaaov", wStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "Psa", wStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "Plv", wStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "Pla", wStore[idx].data(), Nsteps); idx++;


		printMatlabArraySimple( outFile, "t", tStore.data(), Nsteps);

		outFile.close();
	}

}



// print params
void ParamSingleFiberMynardValve::printParams()
{
	printParams( paramOutputFile );
}



void ParamSingleFiberMynardValve::printParams( string outputFileName ){

	if( getProcID() == 0)
	{
		ofstream outFile( outputFileName.c_str() , ios::out );

		// Lump models
		printMatlabVariableSimple( outFile, "Csa", Csa);

		printMatlabVariableSimple( outFile, "Rpv", Rpv);
		printMatlabVariableSimple( outFile, "Rsa", Rsa);

		printMatlabVariableSimple( outFile, "Ppv_fixed", Ppv_fixed);
		printMatlabVariableSimple( outFile, "Psp_fixed", Psp_fixed);
		printMatlabVariableSimple( outFile, "Rpv", Rpv);


		// Single fiber models
		printMatlabVariableSimple( outFile, "V0_lv", V0_lv);
		printMatlabVariableSimple( outFile, "Vw_lv", Vw_lv);
		printMatlabVariableSimple( outFile, "Ta0_lv", Ta0_lv);
		printMatlabVariableSimple( outFile, "Tp0_lv", Tp0_lv);
		printMatlabVariableSimple( outFile, "cp_lv", cp_lv);
		printMatlabVariableSimple( outFile, "cv_lv", cv_lv);
		printMatlabVariableSimple( outFile, "v0_lv", v0_lv);
		printMatlabVariableSimple( outFile, "ls0_lv", ls0_lv);
		printMatlabVariableSimple( outFile, "lsmax_lv", lsmax_lv);
		printMatlabVariableSimple( outFile, "lsw_lv", lsw_lv);


		printMatlabVariableSimple( outFile, "V0_la", V0_la);
		printMatlabVariableSimple( outFile, "Vw_la", Vw_la);
		printMatlabVariableSimple( outFile, "Ta0_la", Ta0_la);
		printMatlabVariableSimple( outFile, "Tp0_la", Tp0_la);
		printMatlabVariableSimple( outFile, "cp_la", cp_la);
		printMatlabVariableSimple( outFile, "cv_la", cv_la);
		printMatlabVariableSimple( outFile, "v0_la", v0_la);
		printMatlabVariableSimple( outFile, "ls0_la", ls0_la);
		printMatlabVariableSimple( outFile, "lsmax_la", lsmax_la);
		printMatlabVariableSimple( outFile, "lsw_la", lsw_la);


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


		// activation model
		printMatlabVariableSimple( outFile, "m1_ventricle", m1_ventricle);
		printMatlabVariableSimple( outFile, "m2_ventricle", m2_ventricle);
		printMatlabVariableSimple( outFile, "tau1_ventricle", tau1_ventricle);
		printMatlabVariableSimple( outFile, "tau2_ventricle", tau2_ventricle);
		printMatlabVariableSimple( outFile, "Ts_ventricle", Ts_ventricle);

		printMatlabVariableSimple( outFile, "m1_atrium", m1_atrium);
		printMatlabVariableSimple( outFile, "m2_atrium", m2_atrium);
		printMatlabVariableSimple( outFile, "tau1_atrium", tau1_atrium);
		printMatlabVariableSimple( outFile, "tau2_atrium", tau2_atrium);
		printMatlabVariableSimple( outFile, "Ts_atrium", Ts_atrium);

		// Other
		printMatlabVariableSimple( outFile, "Tc", Tc);
		printMatlabVariableSimple( outFile, "Ncyc", Ncyc);
		printMatlabVariableSimple( outFile, "Nvars", Nvars);
		printMatlabVariableSimple( outFile, "dt", dt);

		printMatlabVariableSimple( outFile, "xMinDiff", xMinDiff);
		printMatlabVariableSimple( outFile, "fNewtMin", fNewtMin);
		printMatlabVariableSimple( outFile, "dX", dX);
		printMatlabVariableSimple( outFile, "maxNewtIter", maxNewtIter);

		outFile.close();
	}
}


