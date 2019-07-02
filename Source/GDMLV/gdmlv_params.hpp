/*
 * gdmlv_params.hpp
 *
 *  Created on: Sep 11, 2018
 *      Author: brian
 */

#ifndef GDMLV_GDMLV_PARAMS_HPP_
#define GDMLV_GDMLV_PARAMS_HPP_




#include <iostream>
#include <string>
#include <vector>

// Utility library headers
#include "ProlateSplines/prolateSplines.hpp"
#include "CubicSplines1D/cubicSplinesFunctions1D.hpp"

#include "../LVFunctions/LVBicubicDeformation.hpp"
#include "../LVFunctions/LVFourierDeformation.hpp"
#include "../LVFunctions/LVReferenceAdjustments.hpp"

#include "MathematicalOperators/mathematicalOperators.hpp"
#include "ImportExport/importExport.hpp"
#include "DataContainer/dataContainer.hpp"
#include "ParameterParser/ParameterParser.hpp"

using namespace std;


class ParamGdmLV
{

  public:


	int deformationModelIndex;
	int LVSolutionSetting;

	int procID;

	int Nw;

	int N, Nmu, Nnu, Nphi, NmuLid;

	double Tc, Ta;
	double dtMin, dtMax;

	double Ncyc;

	// initial spline shape parameters
	int NnuSplines, NphiSplines;
	int NtopSplines;
	int NSplines;
	int NSplineParams;
	int Ncin0;
	int Ncout0;

	double muin0_const, muout0_const;
	double muSplineConst;

	int Ne0;
	double * e0;
	double * cin0, * cout0;

	double * theta0, * shift0;
	double nuUp0Min;

	double a_original, a;

	// spline deformation parameters
	int muinNnuGrid, muinNphiGrid, nuNnuGrid, nuNphiGrid, phiNnuGrid, phiNphiGrid;

	// deformation parameters
	int Nc, Nd, Ne, Nq;


	double mue;

	double psi_in_b0, psi_out_b0, nuTrans;
	double deltaq;

	double ke, bff, bxx, bfx;
	double kv;
	double ka, kav, ked, kd;
	double Ls0, Lsmax, Lsw, twoLswSq;

	double Rmvo, Rmvc, Raovo, Raovc, beta;
	double Rpao;

	int newtMaxIter;
	double newtXConvg;
	double newtHConvg;

	// Lump parameters
	double Psa_fixed, Ppv_fixed;
	double Cla, Cpv, Cra, Csa, Csp, Csv, Cpa, Cpp;
	double Vu_sa, Vu_sp, Vu_sv, Vu_pa, Vu_pp, Vu_pv, Vu_ra, Vu_rv, Vu_la;
	double Rsa, Rsp, Rsv, Rpa, Rpp, Rpv, Rra;
	double P0_rv, ke_rv, Emax_rv, kr_rv;
	double Vt0;

	double activePower;

	int plotOutput, useRVPressure;
	int NrvBndCoef;

	int rivlinMooneyMaterial;

	// initial spline surfaces contained in classes
	int useRefAdj; // this indcates if the reference adjustment should be used
	// these determine how many modes should be used in the reference shape adjustment:
	int refAdjust_muinNnuBasis, refAdjust_muinNphiBasis;
	LVReferenceShape lvRef; // this replaces the endo0 and epi0
	// ProlateSplines endo0;
	// ProlateSplines epi0;

	// top 1D spline surface class
	Spline1D nuUp0Spline;


	string matOut;
	string initialValuesFile;
	string plotCyclesFile;
	string cartesianOutFile;
	string matlabParamFile;
	string optimParamsFile;
	string optimOutputPrefix;
	string optimOutputPostfix;
	string staticOutputFileName;

	int updateInitialValues;

	// initial values
	double *q0, *w0, *dwdt0, *dqdt0, Plv0, Ppao0;
	double * rvBoundaryCoef;


	// static parameters
	string staticFile;

	double Plv_static, Prv_static, At_static;


	ParamGdmLV(){
		rivlinMooneyMaterial = 0;
		q0 = NULL;
		w0 = NULL;
		dwdt0 = NULL;
		dqdt0 = NULL;
		e0 = NULL;
		rvBoundaryCoef = NULL;
		cin0 = NULL;
		cout0 = NULL;
		theta0 = NULL;
		shift0 = NULL;

	};
	~ParamGdmLV(){
		if(q0 != NULL){ delete [] q0; }
		if(w0 != NULL){ delete [] w0; }
		if(dwdt0 != NULL){ delete [] dwdt0; }
		if(dqdt0 != NULL){ delete [] dqdt0; }
		if(e0 != NULL){ delete [] e0; }
		if(rvBoundaryCoef != NULL){ delete [] rvBoundaryCoef; }
		if(cin0 != NULL){ delete [] cin0; }
		if( cout0 != NULL){ delete [] cout0; }
		if( theta0  != NULL){ delete [] theta0 ; }
		if( shift0  != NULL){ delete [] shift0 ; }

	};

	void computeTopSpline();

	void importParams( DataContainer & prmData );


	void printMatlabParamFile();

};


void computeTopSplineFromVector( double * e0, int NtopSplines, Spline1D * nuUp0Spline, double & nuUp0Min );


#endif /* GDMLV_GDMLV_PARAMS_HPP_ */
