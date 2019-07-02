/*
 * LVReferenceAdjustments.hpp
 *
 *  Created on: May 16, 2019
 *      Author: brian
 */

#ifndef LVFUNCTIONS_LVREFERENCEADJUSTMENTS_HPP_
#define LVFUNCTIONS_LVREFERENCEADJUSTMENTS_HPP_


// library functions
#include "UtilityFunctions/utilityFunctions.hpp"
#include "LVFourierDeformation.hpp"
#include "LVBicubicDeformation.hpp"
#include "LVInitial.hpp"

// namespace
using namespace std;

// local class
class LVReferenceShape {

  protected:
	// Switch that determines if the reference frame adjustments should be used
	bool useRefAdj;

  public:

	// Container for the focal length of the ellipse
	double a_original, delta_a;

	// Containers for the regular surfaces
	ProlateSplines endo0;
	ProlateSplines epi0;

	// Container for the adjustment of the initial surfaces
	FourierDeformation fourierDef0;

	// Containers for the parameters that determine the surfaces
	vector <double> cin0;
	vector <double> cout0;
	vector <double> q0s;

	// size of fourier deformations
	int refAdjust_muinNnuBasis, refAdjust_muinNphiBasis;

	// other locally stored values
	double nuUp0Min;

	LVReferenceShape()
	{
		// The normal setting is to NOT use the reference adjustment.
		// This must be turned on manually by an external function
		useRefAdj = false;
	}
	~LVReferenceShape(){};

	// Class functions
	void refAdjOn(){ useRefAdj = true; }
	void refAdjOff(){ useRefAdj = false; }
	void setupRegularRefSurfaces( vector<double> & cin0, vector <double> & cout0,
			int NnuSplines, int NphiSplines, double nuUp0Min, double muSplineConst, double aIn );
	void setupRefSurfAdj( int refAdjust_muinNnuBasisIn, int refAdjust_muinNphiBasisIn );
	int getRefSurfAdjSize(){ return fourierDef0.Nc + 1; };
	void setRefSurfAdjValues( vector<double> & q0sIn );

	// essential function that will be called (these functions assume that the surfaces have already been setup)
	double focalLength();
	void computeRefSurfaceValuesAndDerivatives( double nu0, double phi0,
			double & muin0, double & muout0, double & dmuin0dnu0, double & dmuin0dphi0, double & dmuout0dnu0, double & dmuout0dphi0 );
	void computeRefSurfaceValuesAdjustment( double nu0, double phi0,
			double muin0, double muout0, double dmuin0dnu0, double dmuin0dphi0, double dmuout0dnu0, double  dmuout0dphi0,
			double & muin0s, double & muout0s, double & dmuin0sdnu0, double & dmuin0sdphi0, double & dmuout0sdnu0, double & dmuout0sdphi0 );
	void computeRefSurfaceValuesOnly( double nu0, double phi0, double & muin0, double & muout0  );
};

// local functions
void hFuncAndDerivatives_refFuncs( double mu0, double nu0, double dmu0dnu0, double dmu0dphi0,
		double & h, double & dhdnu0, double & dhdphi0 );

#endif /* LVFUNCTIONS_LVREFERENCEADJUSTMENTS_HPP_ */
