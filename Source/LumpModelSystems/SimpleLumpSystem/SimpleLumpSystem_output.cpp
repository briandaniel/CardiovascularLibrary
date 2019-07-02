/*
 * simpleLumpSystem_output.cpp
 *
 *  Created on: Sep 14, 2018
 *      Author: brian
 */


#include "SimpleLumpSystem.hpp"







void outputModelResultSimpleLump( int Nsteps, ParamSimpleLump & lumpPrm, double * tStore, double * PlaStore, double * PpvStore,
		double * PppStore, double * PpaStore, double * PraStore, double * PspStore, double * PsaStore, double * VlvStore,
		double * VrvStore, double * PlvStore, double * PrvStore, double * Pmax_lvStore, double * Pmax_rvStore, double * PsvStore,
		double * FilvStore, double * FolvStore, double * FirvStore, double * ForvStore, double * At_lvStore, double * At_rvStore )
{
	ofstream outFile( lumpPrm.outputFile.c_str() , ios::out );

	printMatlabArraySimple( outFile, "t", tStore, Nsteps);
	printMatlabArraySimple( outFile, "Pla", PlaStore, Nsteps);
	printMatlabArraySimple( outFile, "Ppv", PpvStore, Nsteps);
	printMatlabArraySimple( outFile, "Ppv", PpvStore, Nsteps);
	printMatlabArraySimple( outFile, "Ppp", PppStore, Nsteps);
	printMatlabArraySimple( outFile, "Ppa", PpaStore, Nsteps);
	printMatlabArraySimple( outFile, "Pra", PraStore, Nsteps);
	printMatlabArraySimple( outFile, "Psp", PspStore, Nsteps);
	printMatlabArraySimple( outFile, "Psa", PsaStore, Nsteps);
	printMatlabArraySimple( outFile, "Vlv", VlvStore, Nsteps);
	printMatlabArraySimple( outFile, "Vrv", VrvStore, Nsteps);
	printMatlabArraySimple( outFile, "Plv", PlvStore, Nsteps);
	printMatlabArraySimple( outFile, "Prv", PrvStore, Nsteps);
	printMatlabArraySimple( outFile, "Pmax_lv", Pmax_lvStore, Nsteps);
	printMatlabArraySimple( outFile, "Pmax_rv", Pmax_rvStore, Nsteps);
	printMatlabArraySimple( outFile, "Psv", PsvStore, Nsteps);
	printMatlabArraySimple( outFile, "Filv", FilvStore, Nsteps);
	printMatlabArraySimple( outFile, "Folv", FolvStore, Nsteps);
	printMatlabArraySimple( outFile, "Firv", FirvStore, Nsteps);
	printMatlabArraySimple( outFile, "Forv", ForvStore, Nsteps);
	printMatlabArraySimple( outFile, "At_lv", At_lvStore, Nsteps);
	printMatlabArraySimple( outFile, "At_rv", At_rvStore, Nsteps);

	outFile.close();

}



// print params
void ParamSimpleLump::printParams(){

	ofstream outFile( paramOutputFile.c_str() , ios::out );


	printMatlabVariableSimple( outFile, "Cla", Cla );
	printMatlabVariableSimple( outFile, "Cpv", Cpv);
	printMatlabVariableSimple( outFile, "Cpp", Cpp);
	printMatlabVariableSimple( outFile, "Cpa", Cpa);
	printMatlabVariableSimple( outFile, "Cra", Cra);
	printMatlabVariableSimple( outFile, "Csv", Csv);
	printMatlabVariableSimple( outFile, "Csp", Csp);
	printMatlabVariableSimple( outFile, "Csa", Csa);
	printMatlabVariableSimple( outFile, "Vu_la", Vu_la);
	printMatlabVariableSimple( outFile, "Vu_pv", Vu_pv);
	printMatlabVariableSimple( outFile, "Vu_pp", Vu_pp);
	printMatlabVariableSimple( outFile, "Vu_pa", Vu_pa);
	printMatlabVariableSimple( outFile, "Vu_ra", Vu_ra);
	printMatlabVariableSimple( outFile, "Vu_sv", Vu_sv);
	printMatlabVariableSimple( outFile, "Vu_sp", Vu_sp);
	printMatlabVariableSimple( outFile, "Vu_sa", Vu_sa);
	printMatlabVariableSimple( outFile, "Rla", Rla);
	printMatlabVariableSimple( outFile, "Rpv", Rpv);
	printMatlabVariableSimple( outFile, "Rpp", Rpp);
	printMatlabVariableSimple( outFile, "Rpa", Rpa);
	printMatlabVariableSimple( outFile, "Rra", Rra);
	printMatlabVariableSimple( outFile, "Rsv", Rsv);
	printMatlabVariableSimple( outFile, "Rsp", Rsp);
	printMatlabVariableSimple( outFile, "Rsa", Rsa);
	printMatlabVariableSimple( outFile, "P0_rv", P0_rv );
	printMatlabVariableSimple( outFile, "ke_rv", ke_rv);
	printMatlabVariableSimple( outFile, "Emax_rv", Emax_rv);
	printMatlabVariableSimple( outFile, "kr_rv", kr_rv);
	printMatlabVariableSimple( outFile, "Vu_rv", Vu_rv);
	printMatlabVariableSimple( outFile, "P0_lv", P0_lv);
	printMatlabVariableSimple( outFile, "ke_lv", ke_lv);
	printMatlabVariableSimple( outFile, "Emax_lv", Emax_lv);
	printMatlabVariableSimple( outFile, "kr_lv", kr_lv);
	printMatlabVariableSimple( outFile, "Vu_lv", Vu_lv);
	printMatlabVariableSimple( outFile, "Vtotal", Vtotal);
	printMatlabVariableSimple( outFile, "Tc", Tc);
	printMatlabVariableSimple( outFile, "Ta", Ta);
	printMatlabVariableSimple( outFile, "Ncyc", Ncyc);
	printMatlabVariableSimple( outFile, "Nvars", Nvars);
	printMatlabVariableSimple( outFile, "dt", dt);


	outFile.close();
}
