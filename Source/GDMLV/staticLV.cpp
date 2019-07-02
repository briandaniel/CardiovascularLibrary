/*
 * staticLV.cpp
 *
 *  Created on: Dec 7, 2017
 *      Author: brian
 */



#include "gdmlv.hpp"




// assumes that the model has already been evaluated
void GdmLV::computeStaticF( double * F, double Plv ){

	for(int i = 0; i < prm.Nq; i++)
	{
		F[i] = kappa[i] - eta[i]*Plv;
	}

}


void GdmLV::computeStaticLV( double * q, double Plv, double At, double Prv ){

	double t = 0;

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
		localNodes[k].computeStaticStressesTractions( &prm, At );


	}


	//-------------------------------------------------------------------//
	// 4. compute integrals
	//-------------------------------------------------------------------//

	// compute lid integrals (ROOT only)
	lid.computeLidSurfaceIntegrals( etal );
	lid.exportVolumes( Virreg, dVirregdq );




	//-------------------------------------------------------------------//
	// compute integrals
	//-------------------------------------------------------------------//
	computeStaticProcessorIntegrals( localNodes, Nnodes, local_kappa, local_etae,local_Vreg, &prm, Prv );

	// sum the processor integrals
	MPI_Allreduce(local_kappa, kappa, prm.Nq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(local_etae, etae, prm.Nq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&local_Vreg, &Vreg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);



	// need to sum these
	for(int kq = 0; kq < prm.Nq; kq++)
	{
		eta[kq] = etae[kq] + etal[kq];
	}

	// compute the LV volume
	Vlv = Virreg + Vreg;



	// broadcast LV volume
	MPI_Bcast(&Vlv, 1, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(eta, prm.Nq, MPI_DOUBLE, ROOT_ID, MPI_COMM_WORLD);


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


}
























