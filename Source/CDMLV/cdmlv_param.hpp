/*
 * cdmlv_param.hpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */

#ifndef CDMLV_CDMLV_PARAM_HPP_
#define CDMLV_CDMLV_PARAM_HPP_


#include <vector>

#include "ProlateSplines/prolateSplines.hpp"
#include "CubicSplines1D/cubicSplinesFunctions1D.hpp"
#include "../LVFunctions/LVFourierDeformation.hpp"
#include "ImportExport/importExport.hpp"
#include "ParameterParser/ParameterParser.hpp"

class ParamCDMLV{

  public:

	// Domain size
	int Nmu, Nnu, Nphi, NmuLid;

	// Prolate coordinate system size
	double a;

	// Surface variables
	double nuUp0Min;
	int NnuSplines, NphiSplines, NtopSplines, Ne0;
	double muin0_const, muout0_const, muSplineConst;

	// Surface value vectors
	vector<double> e0;
	vector<double> cin0;
	vector<double> cout0;

	// reference frame boundary splines
	ProlateSplines endo0;
	ProlateSplines epi0;
	Spline1D nuUp0Spline;

	// Fiber direction parameters
	double nuTrans;
	double psi_out_b0, psi_in_b0;

	// Deformation parameters
	// Deformable model size
	int Nq;
	int Nc, Nd, Ne;

	// Fourier deformation size
	int muinNnuGrid, muinNphiGrid, nuNnuGrid, nuNphiGrid, phiNnuGrid, phiNphiGrid;
	double dq;

	// Container that stores the number of fourier deformation parameters used
	FourierDeformation fourierDef;

	// Stress parameters
	double ke, bff, bxx, bfx;
	double kv;
	double ka, kav;
	double Ls0, Lsmax, Lsw;

	// Other input files
	string surfaceInputFileName;

	// Imports params
	void importParams( string prmPrefix, DataContainer & prmData );

	// Setup the deformation container and the surface splines
	void setupComputations();

	// Default constructor uses standard values
	ParamCDMLV();

	// More typical constructor requires the prm prefix and container of the parameters
	ParamCDMLV( string prmPrefix, DataContainer & prmData ){
		importParams( prmPrefix, prmData );
	}
};




// other functions
void computeTopSplineFromVector( vector<double> e0, int NtopSplines, Spline1D & nuUp0Spline, double & nuUp0Min );





#endif /* CDMLV_CDMLV_PARAM_HPP_ */







