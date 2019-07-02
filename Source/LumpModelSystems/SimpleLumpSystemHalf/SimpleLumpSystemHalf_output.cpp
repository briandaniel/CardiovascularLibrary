/*
 * Half_output.cpp
 *
 *  Created on: Sep 14, 2018
 *      Author: brian
 */


#include "SimpleLumpSystemHalf.hpp"




void outputModelResultSimpleLumpHalf(int Nsteps, vector<vector<double>> & wStore, vector<vector<double>> & auxStore, vector<double>& tStore, ParamSimpleLumpHalf & lumpPrm )
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
void ParamSimpleLumpHalf::printParams()
{
	printParams( paramOutputFile );
}



void ParamSimpleLumpHalf::printParams( string outputFileName ){

	if( getProcID() == 0)
	{
		ofstream outFile( outputFileName .c_str() , ios::out );

		// Lump models
		printMatlabVariableSimple( outFile, "Csa", Csa);

		printMatlabVariableSimple( outFile, "Rpv", Rpv);
		printMatlabVariableSimple( outFile, "Rsa", Rsa);

		printMatlabVariableSimple( outFile, "Ppv_fixed", Ppv_fixed);
		printMatlabVariableSimple( outFile, "Psp_fixed", Psp_fixed);
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

