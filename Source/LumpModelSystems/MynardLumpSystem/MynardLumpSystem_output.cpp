/*
 * MynardLumpSystem_output.cpp
 *
 *  Created on: Sep 14, 2018
 *      Author: brian
 */


#include "MynardLumpSystem.hpp"







void outputModelResultMynard0D( int Nsteps, ParamMynard0D & lumpPrm, double * tStore, double ** wStore, double ** auxStore, double * wEnd )
{
	ofstream outFile( lumpPrm.outputFile.c_str() , ios::out );

	printMatlabArraySimple( outFile, "wEnd", wEnd, lumpPrm.Nvars);

	// output time array
	printMatlabArraySimple( outFile, "t", tStore, Nsteps);

	// output variables
    // W = [ Vlv, Vrv, Vla, Vra, qmv, qaov, qtri, qpul, zetamv, zetaaov, zetatri, zetapul, Ppv, Ppp, Ppa, Psp, Psa ]
	printMatlabArraySimple( outFile, "Vlv", wStore[0], Nsteps);
	printMatlabArraySimple( outFile, "Vrv", wStore[1], Nsteps);
	printMatlabArraySimple( outFile, "Vla", wStore[2], Nsteps);
	printMatlabArraySimple( outFile, "Vra", wStore[3], Nsteps);
	printMatlabArraySimple( outFile, "qmv", wStore[4], Nsteps);
	printMatlabArraySimple( outFile, "qaov", wStore[5], Nsteps);
	printMatlabArraySimple( outFile, "qtri", wStore[6], Nsteps);
	printMatlabArraySimple( outFile, "qpul", wStore[7], Nsteps);
	printMatlabArraySimple( outFile, "zetamv", wStore[8], Nsteps);
	printMatlabArraySimple( outFile, "zetaaov", wStore[9], Nsteps);
	printMatlabArraySimple( outFile, "zetatri", wStore[10], Nsteps);
	printMatlabArraySimple( outFile, "zetapul", wStore[11], Nsteps);
	printMatlabArraySimple( outFile, "Ppv", wStore[12], Nsteps);
	printMatlabArraySimple( outFile, "Ppp", wStore[13], Nsteps);
	printMatlabArraySimple( outFile, "Ppa", wStore[14], Nsteps);
	printMatlabArraySimple( outFile, "Psp", wStore[15], Nsteps);
	printMatlabArraySimple( outFile, "Psa", wStore[16], Nsteps);

	// output auxiliary variables
	// auxVars = [Plv, Prv, Pla, Pra, Psv, At_ventricle, At_atrium]
	printMatlabArraySimple( outFile, "Plv", auxStore[0], Nsteps);
	printMatlabArraySimple( outFile, "Prv", auxStore[1], Nsteps);
	printMatlabArraySimple( outFile, "Pla", auxStore[2], Nsteps);
	printMatlabArraySimple( outFile, "Pra", auxStore[3], Nsteps);
	printMatlabArraySimple( outFile, "Psv", auxStore[4], Nsteps);
	printMatlabArraySimple( outFile, "At_ventricle", auxStore[5], Nsteps);
	printMatlabArraySimple( outFile, "At_atrium", auxStore[6], Nsteps);


	outFile.close();

}



// print params
void ParamMynard0D::printParams(){

	ofstream outFile( paramOutputFile.c_str() , ios::out );


	// Lump models
	printMatlabVariableSimple( outFile, "Cpv", Cpv);
	printMatlabVariableSimple( outFile, "Cpp", Cpp);
	printMatlabVariableSimple( outFile, "Cpa", Cpa);
	printMatlabVariableSimple( outFile, "Csv", Csv);
	printMatlabVariableSimple( outFile, "Csp", Csp);
	printMatlabVariableSimple( outFile, "Csa", Csa);
	printMatlabVariableSimple( outFile, "Vu_pv", Vu_pv);
	printMatlabVariableSimple( outFile, "Vu_pp", Vu_pp);
	printMatlabVariableSimple( outFile, "Vu_pa", Vu_pa);
	printMatlabVariableSimple( outFile, "Vu_sv", Vu_sv);
	printMatlabVariableSimple( outFile, "Vu_sp", Vu_sp);
	printMatlabVariableSimple( outFile, "Vu_sa", Vu_sa);
	printMatlabVariableSimple( outFile, "Rpv", Rpv);
	printMatlabVariableSimple( outFile, "Rpp", Rpp);
	printMatlabVariableSimple( outFile, "Rpa", Rpa);
	printMatlabVariableSimple( outFile, "Rsv", Rsv);
	printMatlabVariableSimple( outFile, "Rsp", Rsp);
	printMatlabVariableSimple( outFile, "Rsa", Rsa);


	// Varying elastance models
	printMatlabVariableSimple( outFile, "Emax_lv", Emax_lv);
	printMatlabVariableSimple( outFile, "Emin_lv", Emin_lv);
	printMatlabVariableSimple( outFile, "Ks_lv", Ks_lv);
	printMatlabVariableSimple( outFile, "V0_lv", V0_lv);

	printMatlabVariableSimple( outFile, "Emax_rv", Emax_rv);
	printMatlabVariableSimple( outFile, "Emin_rv", Emin_rv);
	printMatlabVariableSimple( outFile, "Ks_rv", Ks_rv);
	printMatlabVariableSimple( outFile, "V0_rv", V0_rv);

	printMatlabVariableSimple( outFile, "Emax_la", Emax_la);
	printMatlabVariableSimple( outFile, "Emin_la", Emin_la);
	printMatlabVariableSimple( outFile, "Ks_la", Ks_la);
	printMatlabVariableSimple( outFile, "V0_la", V0_la);

	printMatlabVariableSimple( outFile, "Emax_ra", Emax_ra);
	printMatlabVariableSimple( outFile, "Emin_ra", Emin_ra);
	printMatlabVariableSimple( outFile, "Ks_ra", Ks_ra);
	printMatlabVariableSimple( outFile, "V0_ra", V0_ra);


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

	printMatlabVariableSimple( outFile, "Aeff_max_tri", Aeff_max_tri);
	printMatlabVariableSimple( outFile, "Aeff_min_tri", Aeff_min_tri);
	printMatlabVariableSimple( outFile, "leff_tri", leff_tri);
	printMatlabVariableSimple( outFile, "Kvo_tri", Kvo_tri);
	printMatlabVariableSimple( outFile, "Kvc_tri", Kvc_tri);
	printMatlabVariableSimple( outFile, "deltaP_open_tri", deltaP_open_tri);
	printMatlabVariableSimple( outFile, "deltaP_close_tri", deltaP_close_tri);

	printMatlabVariableSimple( outFile, "Aeff_max_pul", Aeff_max_pul);
	printMatlabVariableSimple( outFile, "Aeff_min_pul", Aeff_min_pul);
	printMatlabVariableSimple( outFile, "leff_pul", leff_pul);
	printMatlabVariableSimple( outFile, "Kvo_pul", Kvo_pul);
	printMatlabVariableSimple( outFile, "Kvc_pul", Kvc_pul);
	printMatlabVariableSimple( outFile, "deltaP_open_pul", deltaP_open_pul);
	printMatlabVariableSimple( outFile, "deltaP_close_pul", deltaP_close_pul);


	// Other
	printMatlabVariableSimple( outFile, "Vtotal", Vtotal);
	printMatlabVariableSimple( outFile, "Tc", Tc);
	printMatlabVariableSimple( outFile, "Ncyc", Ncyc);
	printMatlabVariableSimple( outFile, "Nvars", Nvars);
	printMatlabVariableSimple( outFile, "dt", dt);


	outFile.close();
}
