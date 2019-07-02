/*
 * cdmlv.hpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */

#ifndef CDMLV_CDMLV_HPP_
#define CDMLV_CDMLV_HPP_


#include "cdmlv_param.hpp"
#include "cdmlv_node.hpp"
#include "cdmlv_integrals.hpp"
#include "cdmlv_lid.hpp"

class CDMLVModel{

  public:

	LidCDMLV lid;
	ParamCDMLV prm;
	vector<CDMNode> localNodes;

	// ints
	int Nnodes;

	// doubles
	double Vlv, Vreg, Virreg;
	double Vlv0;

	// int arrays
	vector<int> local_to_global_id;
	vector<int> proc_node_start;
	vector<int> Nnodes_per_proc;

	// double arrays
	vector<double> dVirregdq;
	vector<double> etal;
	vector<double> kappa;
	vector<double> etae;
	vector<double> dVregdq;
	vector<double> dVlvdq;
	vector<double> eta;
	vector<double> qFinal;
	vector< vector<double> > alpha;


	// setup model
	void setupCDMLVModel( string prmPrefix, DataContainer & prmData );

	// compute dynamic model
	void computeLV( vector<double> & q, double t, double At );

	// static model
	void computeStaticLV( vector<double> & q, double At );
	void computeStaticF( vector<double> & F, double Plv );

	double computeZeroPressureVolume();

	// main constructor
	CDMLVModel( string prmPrefix, DataContainer & prmData )
	{
		setupCDMLVModel( prmPrefix, prmData );
	}

};


#endif /* CDMLV_CDMLV_HPP_ */
