/*
 * gdmlv.cpp
 *
 *  Created on: Sep 10, 2018
 *      Author: brian
 */

#include "gdmlv.hpp"





void GdmLV::computeLV( double * q, double t, double At, double Prv ){



	// temporary local variables
	double * local_kappa = new double [prm.Nq];
	double * local_etae = new double [prm.Nq];
	double * local_dVregdq = new double [prm.Nq];
	double local_Vreg;
	double ** local_alpha = new double * [prm.Nq];
	for(int k = 0; k < prm.Nq; k++)
		local_alpha[k] = new double [prm.Nq];
	double * qdqk = new double [prm.Nq];



	double **** uVecLiddq = new double *** [prm.Nq];
	for(int kq = 0; kq < prm.Nq; kq++)
	{
		uVecLiddq[kq] = new double ** [lid.NmuLid];
		for(int i = 0; i < lid.NmuLid; i++ )
		{
			uVecLiddq[kq][i] = new double * [lid.Nphi];
			for(int j = 0; j < lid.Nphi; j++)
			{
				uVecLiddq[kq][i][j] = new double [3];
				for(int k = 0; k < 3; k++)
				{
					uVecLiddq[kq][i][j][k] = 0.0;
				}
			}
		}
	}
	double * Virreg1_dq = new double [prm.Nq];
	double * Virreg2_dq = new double [prm.Nq];



	//-------------------------------------------------------------------//
	// 1. Compute myocardial kinematics at qdqk
	//-------------------------------------------------------------------//
	// First compute the deformations, etc, at the adjusted q
	for(int kq = 0; kq < prm.Nq; kq++)
	{
		// adjust qdqk
		for (int m = 0; m < prm.Nq; m++)
			qdqk[m] = q[m];
		qdqk[kq] = q[kq] + prm.deltaq;

		// cycle through local nodes computing values
		for (int k = 0; k < Nnodes; k++)
		{

			// compute the deformations and the strains at qdqk
			localNodes[k].deformationStrainComputations( qdqk, fourierDef, &prm );


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

			// cout << k << "  " << localNodes[k].E[0][0] << endl;
		}


		// lid computations, only on ROOT
		lid.computeLid( qdqk, fourierDef, prm );

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
		localNodes[k].deformationStrainComputations( q, fourierDef, &prm );

		// The derivatives are finite difference approximations using the differenced values computed above
		// Note that the variables change meanings during this function to save on storage space
		// First compute the deformations, etc, at the adjusted q
		for(int kq = 0; kq < prm.Nq; kq++)
		{

			localNodes[k].dVIdq[kq] = ( localNodes[k].dVIdq[kq] - localNodes[k].VI )/ prm.deltaq;

			for(int i = 0; i < 3; i++)
			{
				localNodes[k].dudq[kq][i] = ( localNodes[k].dudq[kq][i] - localNodes[k].udisp[i] )/ prm.deltaq;


				for(int j = 0; j < 3; j++)
				{
					localNodes[k].dEdq[kq][i][j] = ( localNodes[k].dEdq[kq][i][j] - localNodes[k].E[i][j] )/prm.deltaq;
					localNodes[k].dEfibdq[kq][i][j] = ( localNodes[k].dEfibdq[kq][i][j] - localNodes[k].Efib[i][j] )/prm.deltaq;
				}
			}

			if( localNodes[k].surfNode == 1 && procID == 0)
			{
			// cout << "k = " << k << ", kq = " << kq << ", dudq =  [" << localNodes[k].dudq[kq][0] << " " << localNodes[k].dudq[kq][1] << " " << localNodes[k].dudq[kq][2] << " ]" << endl;
			}


		}

	}


	//-------------------------------------------------------------------//
	// 3. Compute lid finite difference derivatives
	//-------------------------------------------------------------------//

	// lid computations, redundant on all processors except ROOT
	// this computes and stores the FD derivatives in the "lid" class

	// recompute the lid model at q
	lid.computeLid( q, fourierDef, prm );


	// compute finite difference values using the locally stored values from step 1
	for(int kq = 0; kq < prm.Nq; kq++)
	{


		for(int i = 0; i < lid.Nr; i++)
		{
			for(int j = 0; j < lid.Nphi; j++)
			{
				for(int k = 0; k < 3; k ++)
				{
					lid.duLiddq[kq][i][j][k] = (uVecLiddq[kq][i][j][k] - lid.uLid[i][j][k])/prm.deltaq;
				}
			}
		}

		lid.dVirreg1dq[kq] = (Virreg1_dq[kq] - lid.Virreg1)/prm.deltaq;
		lid.dVirreg2dq[kq] = (Virreg2_dq[kq] - lid.Virreg2)/prm.deltaq;

	}



	//-------------------------------------------------------------------//
	// 4. Compute stresses and surface tractions
	//-------------------------------------------------------------------//
	for (int k = 0; k < Nnodes; k++)
	{
		// Compute stresses
		localNodes[k].computeStressesTractions( &prm, t, At);
	}


	//-------------------------------------------------------------------//
	// 4. compute integrals
	//-------------------------------------------------------------------//

	// compute lid integrals (ROOT only)
	lid.computeLidSurfaceIntegrals( etal );
	lid.exportVolumes( Virreg, dVirregdq );



	// compute myocardial integrals
	computeProcessorIntegrals( localNodes, Nnodes, local_alpha, local_kappa, local_etae, local_dVregdq, local_Vreg, &prm, Prv );

	MPI_Reduce( local_kappa, kappa, prm.Nq, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	MPI_Reduce( local_etae, etae, prm.Nq, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	MPI_Reduce( &local_Vreg, &Vreg, 1, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	MPI_Reduce( local_dVregdq, dVregdq, prm.Nq, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	for(int kq = 0; kq < prm.Nq; kq++)
		MPI_Reduce( local_alpha[kq], alpha[kq], prm.Nq, MPI_DOUBLE, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);


	// need to sum these
	for(int kq = 0; kq < prm.Nq; kq++)
	{
		eta[kq] = etae[kq] + etal[kq];
		// eta[kq] = etae[kq];
	    // eta[kq] = 0;
		dVlvdq[kq] = dVregdq[kq] + dVirregdq[kq];
	}


	// compute the LV volume
	Vlv = Virreg + Vreg;


	// Ensure values are same on all processors
	MPI_Bcast(kappa, prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(etae, prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(&Vreg, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(dVregdq, prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(dVlvdq, prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	for(int kq = 0; kq < prm.Nq; kq++)
		MPI_Bcast(alpha[kq], prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(&Vlv, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);

	if(kappa[0] != kappa[0])
	{


		if( procID == ROOT_ID)
		{
			print2DArray(alpha, prm.Nq, prm.Nq, 3);

			print1DArrayLine(kappa, prm.Nq, 3, " kappa" );
			print1DArrayLine(etae, prm.Nq, 3, " etae" );
			print1DArrayLine(dVregdq, prm.Nq, 3, " dVregdq" );
			print1DArrayLine(etal, prm.Nq, 3, " etal" );
			print1DArrayLine(dVirregdq, prm.Nq, 3, " dVirregdq" );


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

	// clean up temporary variables
	for(int kq = 0; kq < prm.Nq; kq++)
		delete [] local_alpha[kq];
	delete [] local_kappa;
	delete [] local_etae;
	delete [] local_dVregdq;
	delete [] local_alpha;

	delete [] qdqk;
	delete [] Virreg1_dq;
	delete [] Virreg2_dq;


	for(int kq = 0; kq < prm.Nq; kq++)
	{
		for(int i = 0; i < lid.NmuLid; i++ )
		{
			for(int j = 0; j < lid.Nphi; j++)
			{
				delete [] uVecLiddq[kq][i][j];
			}
			delete [] uVecLiddq[kq][i];
		}

		delete [] uVecLiddq[kq];
	}
	delete [] uVecLiddq;


	MPI_Barrier(MPI_COMM_WORLD);
}



double GdmLV::computeZeroPressureVolume(){

	// temporary local variables
	double * q = new double [prm.Nq];
	for(int k = 0; k < prm.Nq; k++)
		q[k] = 0.0;

	// compute model at q = 0, which automatically gives the initial volume.

	// cout << "Computing zero pressure volume. " << endl;
	computeLV( q, 0.0, 0.0, 0.0);

	double Vlv0 = Vlv;


	// clean up temporary variables
	delete [] q;

	//cout << "Zero pressure volume = " << Vlv0 << endl;
	return Vlv0;

}



void GdmLV::updateFiberRotation()
{
	// cycle through local nodes computing values
	for (int k = 0; k < Nnodes; k++)
	{
		localNodes[k].computeFiberRotation( &prm );
	}

}


