/*
 * cdmlv_node.hpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */

#ifndef CDMLV_CDMLV_NODE_HPP_
#define CDMLV_CDMLV_NODE_HPP_

#define DIM 3

#ifndef PI
#define PI 3.141592653589793
#endif

#include <vector>
#include <cmath>
#include "cdmlv_param.hpp"
#include "../LVFunctions/LVInitial.hpp"
#include "../LVFunctions/LVFourierDeformation.hpp"
#include "../LVFunctions/LVDeformation.hpp"
#include "../MuscleFiberFunctions/SimpleMuscleFiber.hpp"


class CDMNode{


  public:

	// local copies of variables
	int Nq;
	double a;
	double acube;

	// load balance variables
	int globalID; // Node global ID (same on all processors)
	int localID; // Node local ID (id for the current processor)

	// coordinate variables
	int imu, jnu, kphi; // mu0, nu0, phi0, coordinate indices
	double mu0, nu0, phi0; // node coordinate in the reference frame
	double nuUp0;
	double umu0, muin0, muout0;
	double x0, y0, z0; // cartesian coordinate

	// surface determinations
	bool endoSurfNode; // Indicator: = 1 if this node is on the endocardial surface

	// precomputed variables
	double snu0, cnu0, snu0sq, cnu0sq;
	double shmu0, chmu0, shmu0sq, chmu0sq;
	double gmu0, gnu0, gphi0, gmult0;
	double chmuin0, chmuin0sq, shmuin0;
	double cphi0, sphi0;
	double fc0, dfc0dnu0, dfc0dphi0;
	double dmuin0dnu0, dmuin0dphi0;
	double dmuout0dnu0, dmuout0dphi0;

	// reference coordinate integral weights
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
	double K;
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

	double lambda; // fiber stretch ratio

	// aux values
	double vCrossTopNorm;
	double vCrossSideNorm;

	// vectors
	std::vector <double> udisp;
	std::vector <double> vcrossSide;
	std::vector <double> vcrossTop;

	// static tensors
	std::vector< std::vector<double> > Q;
	std::vector< std::vector<double> > QT;

	// kinematic tensors
	std::vector< std::vector<double> > F;
	std::vector< std::vector<double> > FT;
	std::vector< std::vector<double> > Finv;
	std::vector< std::vector<double> > C;
	std::vector< std::vector<double> > Cinv;

	std::vector< std::vector<double> > E;
	std::vector< std::vector<double> > Efib;
	std::vector< std::vector<double> > estrain;
	std::vector< std::vector<double> > efib;

	// variables that are required for d/dq of the kinematic tensors
	std::vector <double> dVIdq; // tensor with size Nq
	std::vector< std::vector<double> > dudq; // tensor with size [Nq, 3]
	std::vector< std::vector< std::vector<double> > > dEdq; // tensor with size [Nq, 3, 3]
	std::vector< std::vector< std::vector<double> > > dEfibdq; // tensor with size [Nq, 3, 3]

	// stress tensors
	std::vector< std::vector<double> > Se; // tensor with size [ 3, 3]
	std::vector< std::vector<double> > Sa0; // tensor with size [ 3, 3]
	std::vector< std::vector< std::vector<double> > >  Svq; // tensor with size [Nq, 3, 3]
	std::vector< std::vector< std::vector<double> > >  Saq; // tensor with size [Nq, 3, 3]

	// aux vectors
	std::vector <double> nSide;
	std::vector <double> nTop;
	std::vector <double> tractionForcingTop;
	std::vector <double> tractionForcingSide;
	std::vector <double> bodyForcing;
	std::vector <double> divSigma_cart;
	std::vector< std::vector<double> > sigma_cart;
	std::vector< std::vector<double> > sigma_prol;


	// Default initializer requires the number of deformation modes be defined
	CDMNode( int NqIn ){

		// Only variable set in initializer
		Nq = NqIn;


		// Everything else is -1
		a = -1; acube = -1;
		globalID = -1; localID = -1;
		imu = -1; jnu = -1; kphi = -1;
		mu0 = -1; nu0 = -1; phi0 = -1;
		nuUp0 = -1; umu0 = -1; muin0 = -1; muout0 = -1;
		x0 = -1; y0 = -1; z0 = -1;
	    endoSurfNode = -1;
		snu0 = -1; cnu0 = -1; snu0sq = -1; cnu0sq = -1;
		shmu0 = -1; chmu0 = -1; shmu0sq = -1; chmu0sq = -1;
		gmu0 = -1; gnu0 = -1; gphi0 = -1; gmult0 = -1;
		chmuin0 = -1; chmuin0sq = -1; shmuin0 = -1;
		cphi0 = -1; sphi0 = -1;
		fc0 = -1; dfc0dnu0 = -1; dfc0dphi0 = -1;
		dmuin0dnu0 = -1; dmuin0dphi0 = -1;
		dmuout0dnu0 = -1; dmuout0dphi0 = -1;
		IV0weight = -1;
		IA0weight = -1;
		IA0weight_top = -1;
		mu = -1; nu = -1; phi = -1;
		dnudnu0 = -1; dnudphi0 = -1; dphidnu0 = -1; dphidphi0 = -1;
		sinRatio = -1; K = -1;
		fc = -1; dfcdnu0 = -1; dfcdphi0 = -1;
		d2nudnu02 = -1; d2nudphi0dnu0 = -1; d2nudphi02 = -1;
		d2phidphi02 = -1; d2phidphi0dnu0 = -1; d2phidnu02 = -1;
		x = -1; y = -1; z = -1;
		snu = -1; cnu = -1; snusq = -1; cnusq = -1;
		shmu = -1; chmu = -1; shmusq = -1; chmusq = -1;
		dmudmu0 = -1; dmudnu0 = -1; dmudphi0 = -1;
		gmu = -1; gnu = -1; gphi = -1;
		VI = 0;
		dxdnu0 = 0; dydnu0 = 0; dzdnu0 = 0; dxdphi0 = 0; dydphi0 = 0; dzdphi0 = 0;
		lambda = 0;
		vCrossTopNorm = 0;
		vCrossSideNorm = 0;

		// set vectors to the correct size
		udisp.assign(DIM,0);
		vcrossSide.assign(DIM,0);
		vcrossTop.assign(DIM,0);

		//
		Q.assign(DIM, std::vector<double>(DIM,0) );
		QT.assign(DIM, std::vector<double>(DIM,0) );
		F.assign(DIM, std::vector<double>(DIM,0) );
		FT.assign(DIM, std::vector<double>(DIM,0) );
		Finv.assign(DIM, std::vector<double>(DIM,0) );
		C.assign(DIM, std::vector<double>(DIM,0) );
		Cinv.assign(DIM, std::vector<double>(DIM,0) );
		E.assign(DIM, std::vector<double>(DIM,0) );
		Efib.assign(DIM, std::vector<double>(DIM,0) );
		estrain.assign(DIM, std::vector<double>(DIM,0) );
		efib.assign(DIM, std::vector<double>(DIM,0) );

		dVIdq.assign(Nq,0);
		dudq.assign(Nq,std::vector<double>(DIM,0));
		dEdq.assign(Nq, std::vector< std::vector<double> >( DIM, std::vector<double>(DIM,0) ) );
		dEfibdq.assign(Nq, std::vector< std::vector<double> >( DIM, std::vector<double>(DIM,0) ) );

		Se.assign(DIM,std::vector<double>(DIM,0));
		Sa0.assign(DIM,std::vector<double>(DIM,0));
		Svq.assign(Nq, std::vector< std::vector<double> >( DIM, std::vector<double>(DIM,0) ) );
		Saq.assign(Nq, std::vector< std::vector<double> >( DIM, std::vector<double>(DIM,0) ) );

		nSide.assign(DIM,0);
		nTop.assign(DIM,0);
		tractionForcingTop.assign(DIM,0);
		tractionForcingSide.assign(DIM,0);
		bodyForcing.assign(DIM,0);
		divSigma_cart.assign(DIM,0);

		sigma_cart.assign(DIM,std::vector<double>(DIM,0));
		sigma_prol.assign(DIM,std::vector<double>(DIM,0));

	};

	~CDMNode(){ }

	// Functions in cdmlv_node.cpp
	void setupNodeFromPrescribedLocation( double mu0In, double nu0In, double phi0In, ParamCDMLV & prm );
	void setPrescribedLocation( double mu0In, double nu0In, double phi0In, ParamCDMLV & prm );
	void setupNodeFromLoadBalance(int globalID_in, int localID_in, ParamCDMLV & prm );
	void setIndicesFromLoadBalance(int globalID_in, int localID_in, ParamCDMLV & prm );
	void setLocationFromIndices( ParamCDMLV & prm );
	void initialComputations( ParamCDMLV & prm );
	void computeIntegrationWeights( ParamCDMLV & prm );


	// Functions in cdmlv_node_fiber_directions.cpp
	void computeFiberRotation( double * es, double * en, double * ef );
	void computeFiberRotation( ParamCDMLV & prm );
	void computeFiberVectors( double * e1, double * e2, double * e3, double * e4, double * es, double * en, double * ef, ParamCDMLV & prm );

	// Kinematic functions in cdmlv_node_kinematics.cpp
	void computeSubvalues();
	void computeDeformation( vector <double> & q, ParamCDMLV & prm );
	void computeFourierDeformation( vector <double> & q, ParamCDMLV & prm );
	void computeDeformationGradient();
	void computeStrainTensors();
	void rotateProlateToFiber( vector< vector<double> > & T, vector< vector<double> > & Tfib );
	void rotateFiberToProlate( vector< vector<double> > & Tfib, vector< vector<double> > &  T );
	void computeCavityVolume( ParamCDMLV & prm );
	void kinematicComputations( vector <double> & q, ParamCDMLV & prm );

	// Stress functions in cdmlv_node_stress.cpp
	void computeStressesTractions( ParamCDMLV & prm, double t, double At );
	void computeStaticStressesTractions( ParamCDMLV & prm, double At );
	void computeTopSurfaceNorm( ParamCDMLV & prm );
	void surfaceTractions( ParamCDMLV & prm );
	void elasticStress( ParamCDMLV & prm );
	void viscousStress( ParamCDMLV & prm );
	void activeStress0( ParamCDMLV & prm, double At );
	void activeStressTD( ParamCDMLV & prm, double At );
	void computeCauchyStress();

	// Output node values in cdmlv_node_output.cpp
	void printNodeToMFile( ofstream &fileID );

};






#endif /* CDMLV_CDMLV_NODE_HPP_ */
