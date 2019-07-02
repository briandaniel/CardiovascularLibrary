/*
 * gdmlv.hpp
 *
 *  Created on: Sep 10, 2018
 *      Author: brian
 */

#ifndef GDMLV_GDMLV_HPP_
#define GDMLV_GDMLV_HPP_

using namespace std;

// model includes
#include "gdmlv_params.hpp"
#include "gdmlv_node.hpp"
#include "gdmlv_lid.hpp"
#include "equilibriumIntegralFunctions.hpp"

class GdmLV{

  public:


	//----------Local class variables----------//
	// ints
	int Nnodes;
	int procID;
	int Nprocs;
	int maxNodesPerProc;
	// doubles
	double Vreg;
	double Virreg;
	double Vlv;
	double Vlv0;

	// int arrays
	int * local_to_global_id;
	int * proc_node_start;
	int * Nnodes_per_proc;
	// double arrays
	double * dVirregdq;
	double * etal;
	double * kappa;
	double * etae;
	double * dVregdq;
	double * dVlvdq;
	double * eta;
	double * qFinal;
	double ** alpha;

	// class containers
	ParamGdmLV prm;
	FourierDeformation fourierDef;
	vector <NodeGdmLV> localNodes;
	LidGdmLV lid;


	// 1. gdmlv.cpp
	// 1a. compute the ventricle
	void computeLV( double * q, double t, double At, double Prv );
	// 1b. zero pressure volume computation
	double computeZeroPressureVolume();

	// 2. gdmlv_setup_destroy.cpp
	// 2a. sets the lv based on settings in paramLVFile
	void setupLV( DataContainer & prmData  );
	// 2b. In situations where the reference shape has changes
	// the initial locations should be recomputed
	void recomputeInitialLocations();
	// 2c. clears the variables
	void clearLV();

	// 3. staticLV.cpp
	void computeStaticF( double * F, double Plv );
	void computeStaticLV( double * q, double Plv, double At, double Prv );

	// 4. Update the fiber angles (if new psi's were assigned)
	void updateFiberRotation();

};


// other function definitions


#endif /* GDMLV_GDMLV_HPP_ */
