/*
 * gdmlv_setup_destroy.cpp
 *
 *  Created on: Sep 12, 2018
 *      Author: brian
 */



#include "gdmlv.hpp"



void GdmLV::setupLV( DataContainer & prmData ){

	int nodesComputed;

	prm.importParams( prmData );

	//-------------------------------------------------------------------//
	// 1. set up the deformation parameters
	//-------------------------------------------------------------------//
	fourierDef.setSize( prm.muinNnuGrid, prm.muinNphiGrid, prm.nuNnuGrid, prm.nuNphiGrid, prm.phiNnuGrid, prm.phiNphiGrid );

	//-------------------------------------------------------------------//
	// 2. distribute the computation across available processors
	//-------------------------------------------------------------------//

	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	// Distribute load to nodes
	maxNodesPerProc = (int) ceil( (double) prm.N/ (double) Nprocs );
	proc_node_start = new int [Nprocs];

	// Separate load onto the available processors
	Nnodes_per_proc = new int[Nprocs];
	nodesComputed = 0;
	for (int k = 0; k < Nprocs-1; k++)
	{
		Nnodes_per_proc[k] = maxNodesPerProc;
		proc_node_start[k] = nodesComputed;
		nodesComputed = nodesComputed + Nnodes_per_proc[k];
	}
	// Last processor will possibly have slightly fewer nodes
	proc_node_start[Nprocs-1] = nodesComputed;
	Nnodes_per_proc[Nprocs-1] = prm.N-nodesComputed;
	Nnodes = Nnodes_per_proc[procID];
	local_to_global_id = new int [Nnodes];

	for(int k_loc = 0; k_loc < Nnodes_per_proc[procID]; k_loc++)
	{
		local_to_global_id[k_loc] = proc_node_start[procID] + k_loc;
	}

	//-------------------------------------------------------------------//
	// 3. set up myocardial domain
	//-------------------------------------------------------------------//

	// Create a vector of pointers to the local nodes to be computed on this processor
	NodeGdmLV tempNode;

	for(int k = 0; k < Nnodes; k++)
	{
		localNodes.push_back(tempNode);

		// Initialize the nodes and the static values
		localNodes[k].initializeTensors( &prm );
		localNodes[k].setLocation( local_to_global_id[k], k, &prm );
		localNodes[k].initialComputations( &prm );

	}




	//-------------------------------------------------------------------//
	// 4. set up lid
	//-------------------------------------------------------------------//
	// redundant on all processors except root
	// copy all of the nodes onto the main array by summation
	lid.createLid( prm );
	lid.setInitialPositions( prm );




	//-------------------------------------------------------------------//
	// 5. set up variables
	//-------------------------------------------------------------------//
	dVirregdq = new double [prm.Nq];
	etal = new double [prm.Nq];

	// gather integrals to root processor
	kappa = new double [prm.Nq];
	etae = new double [prm.Nq];
	dVregdq = new double [prm.Nq];
	alpha = new double * [prm.Nq];
	for(int k = 0; k < prm.Nq; k++)
		alpha[k] = new double [prm.Nq];
	eta = new double [prm.Nq];
	dVlvdq = new double [prm.Nq];

	qFinal = new double [prm.Nq];


	// compute zero pressure volume
	Vlv0 = computeZeroPressureVolume();


}


// In situations where the reference shape has changes
// the initial locations should be recomputed
void GdmLV::recomputeInitialLocations()
{
	if( prm.useRefAdj == 1)
	{
		prm.a = prm.lvRef.focalLength();
	}

	for(int k = 0; k < Nnodes; k++)
	{
		localNodes[k].setLocation( local_to_global_id[k], k, &prm );
		localNodes[k].initialComputations( &prm );
	}
	lid.setInitialPositions( prm );

	// compute zero pressure volume
	Vlv0 = computeZeroPressureVolume();
}


void GdmLV::clearLV(){



	// clean up
	delete [] proc_node_start;
	delete [] Nnodes_per_proc;
	delete [] local_to_global_id;


	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	// remove the lid
	lid.destroyLid();


	// Delete the nodes (this is required to free the matrix arrays)
	for(int k = 0; k < Nnodes; k++)
	{
		localNodes[k].destroyTensors();
	}

	delete [] dVirregdq;
	delete [] etal;
	for(int k = 0; k < prm.Nq; k++)
		delete [] alpha[k];
	delete [] kappa;
	delete [] etae;
	delete [] dVregdq;
	delete [] alpha;
	delete [] eta;
	delete [] dVlvdq;

	delete [] qFinal;

}

