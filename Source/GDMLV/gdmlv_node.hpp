/*
 * gdmlv_node.hpp
 *
 *  Created on: Sep 11, 2018
 *      Author: brian
 */

#ifndef GDMLV_GDMLV_NODE_HPP_
#define GDMLV_GDMLV_NODE_HPP_


#include <iostream>
#include <string>
#include <vector>

// local headers
#include "gdmlv_params.hpp"
#include "rvPressureFunctions.hpp"

// Utility library headers
#include "ProlateSplines/prolateSplines.hpp"
#include "CubicSplines1D/cubicSplinesFunctions1D.hpp"
#include "../LVFunctions/LVBicubicDeformation.hpp"
#include "../LVFunctions/LVFourierDeformation.hpp"
#include "../LVFunctions/LVDeformation.hpp"
#include "../LVFunctions/LVInitial.hpp"

#include "MathematicalOperators/mathematicalOperators.hpp"
#include "ImportExport/importExport.hpp"
#include "UtilityFunctions/utilityFunctions.hpp"

class NodeGdmLV{

  public:

	// local copies of variables
	double a;
	double acube;
	int Nq;

	// Static NodeGdmLV values
	int globalID; // Global ID
	int localID; // ID on the local processor
	int imu, jnu, kphi; // mu0, nu0, phi0, coordinate indices

	bool surfNode;
	bool endoSurfNode; // Indicator: = 1 if this node is on the endocardial surface
	bool rvSurfNode;

	double mu0, nu0, phi0; // node coordinate in the reference frame
	double nuUp0;
	double umu0, muin0, muout0;

	double chmue;

	double x0, y0, z0; // cartesian coordinate

	double snu0, cnu0, snu0sq, cnu0sq;
	double shmu0, chmu0, shmu0sq, chmu0sq;
	double gmu0, gnu0, gphi0, gmult0;
	double chmuin0, chmuin0sq, shmuin0;
	double cphi0, sphi0;
	double fc0, dfc0dnu0, dfc0dphi0;
	double dmuin0dnu0, dmuin0dphi0;
	double dmuout0dnu0, dmuout0dphi0;

	double K;


	// volume integral weight
	double IV0weight;
	double IA0weight;
	double IA0weight_top;


	// ------------------------------------------------------------------- //
	// Values that depend on the state -- these variables get reused
	// for multiple computations, so all of the computations have to be done
	// in the correct order.
	double mu, nu, phi;
	double dnudnu0, dnudphi0, dphidnu0, dphidphi0;
	double sinRatio;
	double fc, dfcdnu0, dfcdphi0;
	double d2nudnu02, d2nudphi0dnu0, d2nudphi02;
	double d2phidphi02, d2phidphi0dnu0, d2phidnu02;
	double x, y, z;
	double snu, cnu, snusq, cnusq;
	double shmu, chmu, shmusq, chmusq;
	double dmudmu0, dmudnu0, dmudphi0;
	double gmu, gnu, gphi;
	double VI;
	double dxdnu0, dydnu0, dzdnu0, dxdphi0, dydphi0, dzdphi0;

	double eff_ed; // end diastolic stress
	double lambda; // fiber stretch ratio

	// aux values
	double vCrossTopNorm;
	double vCrossSideNorm;

	// vectors
	double * udisp;
	double * vcrossSide;
	double * vcrossTop;

	// static tensors
	double ** Q;
	double ** QT;

	// kinematic tensors
	double ** F;
	double ** FT;
	double ** Finv;
	double ** C;
	double ** Cinv;

	double ** E;
	double ** Efib;
	double ** estrain;
	double ** efib;

	// variables that are required for d/dq of the kinematic tensors
	double * dVIdq; // tensor with size Nq
	double ** dudq; // tensor with size [Nq, 3]
	double *** dEdq; // tensor with size [Nq, 3, 3]
	double *** dEfibdq; // tensor with size [Nq, 3, 3]

	// stress tensors
	double ** Se;
	double ** Sa0; // tensor with size [Nq, 3, 3]
	double *** Svq; // tensor with size [Nq, 3, 3]
	double *** Saq; // tensor with size [Nq, 3, 3]

	// aux vectors
	double * nSide;
	double * nTop;
	double * tractionForcingTop;
	double * tractionForcingSide;
	double * bodyForcing;
	double * divSigma_cart;
	double ** sigma_cart;
	double ** sigma_prol;

	// ------------------------------------------------------------------- //


	NodeGdmLV(){

		// set pointers to null initially
		udisp = NULL;
		vcrossSide = NULL;
		vcrossTop = NULL;
		Q = NULL;
		QT  = NULL;
		F = NULL;
		FT = NULL;
		Finv = NULL;
		C = NULL;
		Cinv = NULL;
		E = NULL;
		Efib = NULL;
		estrain = NULL;
		efib = NULL;
		dVIdq = NULL;
		dudq = NULL;
		dEdq = NULL;
		dEfibdq = NULL;
		Se = NULL;
		Sa0 = NULL;
		Svq = NULL;
		Saq = NULL;
		tractionForcingSide = NULL;
		tractionForcingTop = NULL;
		nSide = NULL;
		nTop = NULL;
		bodyForcing = NULL;
		divSigma_cart = NULL;
		sigma_cart = NULL;
		sigma_prol = NULL;

	};

	void initializeTensors( ParamGdmLV * prm);
	void destroyTensors();



	// node.cpp
	void setLocation(int globalID_in, int localID_in, ParamGdmLV * prm);
	void setPrescribedLocation( double mu0In, double nu0In, double phi0In, int globalID_in, int localID_in, ParamGdmLV * prm);
	void setFixedLocation( int imuIn, int jnuIn, int kphiIn, int Nmu, int Nnu, int Nphi, ParamGdmLV * prm);
	void initialComputations( ParamGdmLV * prm );
	void computeIntegrationWeights( ParamGdmLV * prm );
	void computeSubvalues( );
	void computeSplineDeformation( double * q, ProlateSplines & nuSpline, ProlateSplines & phiSpline, ProlateSplines & fsSpline, ParamGdmLV * prm );
	void printNode( ofstream &fileID  );

	// initialFunctions.cpp
	void computeFiberRotation( double * es, double * en, double * ef ); // using given fiber directions
	void computeFiberRotation( ParamGdmLV * prm ); // default
	void computeFiberVectors( double * e1, double * e2, double * e3, double * e4, double * es, double * en, double * ef, ParamGdmLV * prm );

	// kinematic.cpp
	void computeDeformation( double * q, FourierDeformation & fourierDef, ParamGdmLV * prm  );
	void computeFourierDeformation( double * q, FourierDeformation & fourierDef );
	void computeNodeSevenParamDeformation( double * q );
	void computeDeformationGradient();
	void computeStrainTensors();
	void rotateProlateToFiber( double** T, double** Tfib );
	void rotateFiberToProlate( double** Tfib, double** T );
	void computeDisplacements( );
	void deformationStrainComputations( double * q, FourierDeformation & fourierDef, ParamGdmLV * prm );
	void computeCavityVolume( ParamGdmLV * prm );

	// stress.cpp
	void computeStaticStressesTractions( ParamGdmLV * prm, double At );
	void computeTopSurfaceNorm( ParamGdmLV * prm );
	void surfaceTractions( ParamGdmLV * prm );
	void elasticStress( ParamGdmLV * prm );
	void viscousStress( ParamGdmLV * prm );
	void activeStress0( ParamGdmLV * prm, double At );
	void activeStressTD( ParamGdmLV * prm, double At );
	void computeStressesTractions( ParamGdmLV * prm, double t, double At );
	void computeCauchyStress();

	void elasticStressRivlinMooney( ParamGdmLV * prm );

};



#endif /* GDMLV_GDMLV_NODE_HPP_ */
