/*
 * SingleFiberVarResSystem_output.cpp
 *
 *  Created on: Dec 11, 2018
 *      Author: brian
 */


#include "SingleFiberVarResSystem.hpp"




void outputModelResultSingleFiberVarRes(int Nsteps, vector<vector<double>> & wStore,
		vector<vector<double>> & auxStore, vector<double>& tStore, ParamSingleFiberVarRes & lumpPrm )
{

	if( getProcID() == 0)
	{
		int idx;

		ofstream outFile( lumpPrm.outputFile.c_str() , ios::out );

		// W = [ Vla, Vlv, Psa ]
		idx = 0;
		printMatlabArraySimple( outFile, "Vla", wStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "Vlv", wStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "Psa", wStore[idx].data(), Nsteps); idx++;

		// auxVars = [Plv, Pla, At_ventricle, At_atrium, qmv, qaov, Rmv, Raov ]
		idx = 0;
		printMatlabArraySimple( outFile, "Plv", auxStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "Pla", auxStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "At_ventricle", auxStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "At_atrium", auxStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "qmv", auxStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "qaov", auxStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "Rmv", auxStore[idx].data(), Nsteps); idx++;
		printMatlabArraySimple( outFile, "Raov", auxStore[idx].data(), Nsteps); idx++;

		printMatlabArraySimple( outFile, "t", tStore.data(), Nsteps);

		outFile.close();
	}

}



// print params
void ParamSingleFiberVarRes::printParams()
{
	printParams( paramOutputFile );
}



void ParamSingleFiberVarRes::printParams( string outputFileName ){

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

		// Valve model
		printMatlabVariableSimple( outFile, "Rmvo", Rmvo );
		printMatlabVariableSimple( outFile, "Rmvc", Rmvc);
		printMatlabVariableSimple( outFile, "Raovo", Raovo);
		printMatlabVariableSimple( outFile, "Raovc", Raovc);
		printMatlabVariableSimple( outFile, "beta", beta );
		printMatlabVariableSimple( outFile, "Rpao", Rpao );

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

		printMatlabVariableSimple( outFile, "activeDurationAdjustment", activeDurationAdjustment);


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
