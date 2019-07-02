/*
 * cdmlv.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */


#include "cdmlv.hpp"

void CDMLVModel::setupCDMLVModel( string prmPrefix, DataContainer & prmData ){


	//-------------------------------------------------------------------//
	// 1. load the parameters
	//-------------------------------------------------------------------//
	prm.importParams( prmPrefix, prmData );


	//-------------------------------------------------------------------//
	// 2. distribute the computation across available processors
	//-------------------------------------------------------------------//
	int Nprocs, procID;
	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	// Distribute load to nodes
	int N = prm.Nmu*prm.Nnu*prm.Nphi;
	int maxNodesPerProc = (int) ceil( (double) N/ (double) Nprocs );

	// Separate load onto the available processors
	int nodesComputed = 0;
	for (int k = 0; k < Nprocs-1; k++)
	{
		Nnodes_per_proc.push_back( maxNodesPerProc );
		proc_node_start.push_back( nodesComputed );
		nodesComputed = nodesComputed + Nnodes_per_proc[k];
	}
	// Last processor will possibly have slightly fewer nodes
	proc_node_start.push_back( nodesComputed );
	Nnodes_per_proc.push_back( N - nodesComputed );
	Nnodes = Nnodes_per_proc[procID];

	for(int k_loc = 0; k_loc < Nnodes_per_proc[procID]; k_loc++)
	{
		local_to_global_id.push_back( proc_node_start[procID] + k_loc );
	}

	//-------------------------------------------------------------------//
	// 3. set up myocardial domain
	//-------------------------------------------------------------------//

	// Create a vector of pointers to the local nodes to be computed on this processor
	CDMNode node( prm.Nq );
	localNodes.assign(Nnodes, node );
	for(int k = 0; k < Nnodes; k++)
	{
		// Initialize the nodes and the static values
		localNodes[k].setupNodeFromLoadBalance(local_to_global_id[k], k, prm );
	}


	//-------------------------------------------------------------------//
	// 4. set up lid
	//-------------------------------------------------------------------//
	// redundant on all processors except root
	// copy all of the nodes onto the main array by summation
	lid.setupLid( prm );


	//-------------------------------------------------------------------//
	// 5. set up variables
	//-------------------------------------------------------------------//
	dVirregdq.assign(prm.Nq,0);
	etal.assign(prm.Nq,0);
	kappa.assign(prm.Nq,0);
	etae.assign(prm.Nq,0);
	dVregdq.assign(prm.Nq,0);
	eta.assign(prm.Nq,0);
	dVlvdq.assign(prm.Nq,0);
	qFinal.assign(prm.Nq,0);
	alpha.assign(prm.Nq, std::vector<double>(prm.Nq,0) );


	//-------------------------------------------------------------------//
	// 6. compute zero pressure volume
	//-------------------------------------------------------------------//
	Vlv0 = computeZeroPressureVolume();


}






void CDMLVModel::computeLV( vector<double> & q, double t, double At ){



	// local integral variables
	double local_Vreg;
	vector<double> local_kappa( prm.Nq, 0);
	vector<double> local_etae( prm.Nq, 0);
	vector<double> local_dVregdq( prm.Nq, 0);
	vector< vector<double> > local_alpha( prm.Nq, vector<double>(prm.Nq, 0));
	vector<double> qdqk( prm.Nq, 0);


	// lid definitions
	vector< vector< vector<vector<double>> > > uVecLiddq
		(prm.Nq, vector<vector<vector<double>>>(lid.NmuLid,vector<vector<double>>(lid.Nphi, vector<double>(3,0))));
	vector<double> Virreg1_dq( prm.Nq, 0);
	vector<double> Virreg2_dq( prm.Nq, 0);




	//-------------------------------------------------------------------//
	// 1. Compute myocardial kinematics at qdqk
	//-------------------------------------------------------------------//
	// First compute the deformations, etc, at the adjusted q
	for(int kq = 0; kq < prm.Nq; kq++)
	{
		// adjust qdqk
		for (int m = 0; m < prm.Nq; m++)
			qdqk[m] = q[m];
		qdqk[kq] = q[kq] + prm.dq;

		// cycle through local nodes computing values
		for (int k = 0; k < Nnodes; k++)
		{

			// compute the deformations and the strains at qdqk
			localNodes[k].kinematicComputations( qdqk, prm );

			// Save the strains/positions -- these arent actually the derivatives yet
			localNodes[k].dVIdq[kq] = localNodes[k].VI;

			for(int i = 0; i < 3; i++)
			{
				localNodes[k].dudq[kq][i] = localNodes[k].udisp[i];
				for(int j = 0; j < 3; j++)
				{
					localNodes[k].dEdq[kq][i][j] = localNodes[k].E[i][j];
					localNodes[k].dEfibdq[kq][i][j] = localNodes[k].Efib[i][j];
				}
			}
		}

		// lid computations, only on ROOT
		lid.computeLid( qdqk, prm );

		// compute displacement at dq
		for(int i = 0; i < lid.Nr; i++)
		{
			for(int j = 0; j < lid.Nphi; j++)
			{

				uVecLiddq[kq][i][j][0] = lid.x[i][j] - lid.x0[i][j];
				uVecLiddq[kq][i][j][1] = lid.y[i][j] - lid.y0[i][j];
				uVecLiddq[kq][i][j][2] = lid.z[i][j] - lid.z0[i][j];

			}
		}
		Virreg1_dq[kq] = lid.Virreg1;
		Virreg2_dq[kq] = lid.Virreg2;


	}



	//-------------------------------------------------------------------//
	// 2. Compute myocardial kinematics at q and finite difference derivatives
	//-------------------------------------------------------------------//

	// cycle through local nodes computing values
	for (int k = 0; k < Nnodes; k++)
	{

		// compute the deformations and the strains at q
		localNodes[k].kinematicComputations(q, prm);

		// The derivatives are finite difference approximations using the differenced values computed above
		// Note that the variables change meanings during this function to save on storage space
		// First compute the deformations, etc, at the adjusted q
		for(int kq = 0; kq < prm.Nq; kq++)
		{

			localNodes[k].dVIdq[kq] = ( localNodes[k].dVIdq[kq] - localNodes[k].VI )/ prm.dq;

			for(int i = 0; i < 3; i++)
			{
				localNodes[k].dudq[kq][i] = ( localNodes[k].dudq[kq][i] - localNodes[k].udisp[i] )/ prm.dq;


				for(int j = 0; j < 3; j++)
				{
					localNodes[k].dEdq[kq][i][j] = ( localNodes[k].dEdq[kq][i][j] - localNodes[k].E[i][j] )/prm.dq;
					localNodes[k].dEfibdq[kq][i][j] = ( localNodes[k].dEfibdq[kq][i][j] - localNodes[k].Efib[i][j] )/prm.dq;
				}
			}


		}

	}



	//-------------------------------------------------------------------//
	// 3. Compute lid finite difference derivatives
	//-------------------------------------------------------------------//
	// lid computations, redundant on all processors except ROOT
	// this computes and stores the FD derivatives in the "lid" class

	// recompute the lid model at q
	lid.computeLid( q, prm );


	// compute finite difference values using the locally stored values from step 1
	for(int kq = 0; kq < prm.Nq; kq++)
	{


		for(int i = 0; i < lid.Nr; i++)
		{
			for(int j = 0; j < lid.Nphi; j++)
			{
				for(int k = 0; k < 3; k ++)
				{
					lid.duLiddq[kq][i][j][k] = (uVecLiddq[kq][i][j][k] - lid.uLid[i][j][k])/prm.dq;
				}
			}
		}

		lid.dVirreg1dq[kq] = (Virreg1_dq[kq] - lid.Virreg1)/prm.dq;
		lid.dVirreg2dq[kq] = (Virreg2_dq[kq] - lid.Virreg2)/prm.dq;

	}



	//-------------------------------------------------------------------//
	// 4. Compute stresses and surface tractions
	//-------------------------------------------------------------------//
	for (int k = 0; k < Nnodes; k++)
	{
		// Compute stresses
		localNodes[k].computeStressesTractions( prm, t, At);

	}


	//-------------------------------------------------------------------//
	// 4. compute integrals
	//-------------------------------------------------------------------//

	// compute lid integrals (ROOT only)
	lid.computeLidSurfaceIntegrals( etal );
	lid.exportVolumes( Virreg, dVirregdq );

	// compute regular myocardial integrals
	computeRegularIntegrals( localNodes, Nnodes, local_alpha, local_kappa, local_etae, local_dVregdq, local_Vreg, prm );

	// Reduce all processors to root
	MPI_Reduce( local_kappa.data(), kappa.data(), prm.Nq, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	MPI_Reduce( local_etae.data(), etae.data(), prm.Nq, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	MPI_Reduce( &local_Vreg, &Vreg, 1, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	MPI_Reduce( local_dVregdq.data(), dVregdq.data(), prm.Nq, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	for(int kq = 0; kq < prm.Nq; kq++)
		MPI_Reduce( local_alpha[kq].data(), alpha[kq].data(), prm.Nq, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);


	// need to sum these
	for(int kq = 0; kq < prm.Nq; kq++)
	{
		eta[kq] = etae[kq] + etal[kq];
		dVlvdq[kq] = dVregdq[kq] + dVirregdq[kq];
	}

	// compute the LV volume
	Vlv = Virreg + Vreg;


	// Ensure values are same on all processors
	MPI_Bcast(kappa.data(), prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(etae.data(), prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(eta.data(), prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(&Vreg, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(dVregdq.data(), prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(dVlvdq.data(), prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	for(int kq = 0; kq < prm.Nq; kq++)
		MPI_Bcast(alpha[kq].data(), prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(&Vlv, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);


	if(kappa[0] != kappa[0])
	{


		if( getProcID() == ROOT_ID)
		{
			cout << "CDMLV solver failed..." << endl << endl;

			cout << "   alpha = "; print2DVector(alpha);
			cout << "   kappa = "; print1DVector(kappa);
			cout << "   etae = "; print1DVector(etae);
			cout << "   dVregdq = "; print1DVector(dVregdq);
			cout << "   etal = "; print1DVector(etal);
			cout << "   dVirregdq = "; print1DVector(dVirregdq);

			cout << "-----------------------------------------------------------" << endl << endl;
		}
		// exit(0);

	}

/*
	string fileName = "cCodeOutput.txt";
			outputNodeValues( fileName, local_to_global_id, procID, ROOT_ID, Nnodes, localNodes,  &prm,16 );

	fileName = "cCodeOutput2.txt";
	outputLidVariable( fileName, local_to_global_id, procID,
						ROOT_ID, Nnodes, &lid, &prm, 16 );

	string fileName = "cCodeOutput.txt";
			outputNodeValues( fileName, local_to_global_id, procID, ROOT_ID, Nnodes, localNodes,  &prm,16 );

 */



	MPI_Barrier(MPI_COMM_WORLD);
}



double CDMLVModel::computeZeroPressureVolume(){

	// compute model at q = 0, which automatically gives the initial volume.
	vector<double> q (prm.Nq,0);

	// cout << "Computing zero pressure volume. " << endl;
	computeLV( q, 0.0, 0.0 );

	double Vlv0 = Vlv;

	//cout << "Zero pressure volume = " << Vlv0 << endl;
	return Vlv0;

}
