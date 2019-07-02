/*
 * gdmlv_output.cpp
 *
 *  Created on: Nov 27, 2018
 *      Author: brian
 */

#include "gdmlv_output.hpp"



void printStaticLVShapeMatlabFormat( double * q, string outputName, DataContainer & prmData)
{

	// initialize model
	GdmLV lv;
	lv.setupLV( prmData );


	// visualize model
	int Nmu = 10;
	int Nnu = 20;
	int Nphi = 20;
	int NmuLid = 20;
	int NphiLid = Nphi;


	// Compute the rigid motion rotation matrices
	vector <double> theta(3,0);
	vector <double> shift(3,0);

	string initSurfacePrefix = "initSurf";
	readPrmVector("theta",  initSurfacePrefix, prmData, theta);
	readPrmVector("shift",  initSurfacePrefix, prmData, shift);

	RigidMotion rm;
	rm.computeRigidPrm( theta.data(), shift.data() );


	printModelStateMatlabFormat( q, Nmu, Nnu, Nphi, NmuLid, NphiLid, lv, rm, outputName );

	// clean up
	lv.clearLV();

}





void printLVShapesTimeSeries( vector<vector<double>> & qStore, int Nsteps, int Nskip, string outputFolder, DataContainer & prmData)
{

	// initialize model
	GdmLV lv;
	lv.setupLV( prmData );

	// visualize model
	int Nmu = 7;
	int Nnu = 15;
	int Nphi = 15;
	int NmuLid = 15;
	int NphiLid = Nphi;

	// Compute the rigid motion rotation matrices
	vector <double> theta(3,0);
	vector <double> shift(3,0);
	vector<double> q(lv.prm.Nq,0);

	RigidMotion rm;
	rm.computeRigidPrm( theta.data(), shift.data() );

	for(int i = 0; i < Nsteps; i = i + 1 + Nskip )
	{


		string outputName = outputFolder + "/img" + to_string(i+1) + ".m";
		string outputName2 = outputFolder + "/vals" + to_string(i+1) + ".m";

		if(getProcID() == 0 )
			cout << "Printing surfaces to " << outputName << endl;

		for(int j = 0; j < lv.prm.Nq; j++)
		{
			q[j] = qStore[j][i];
		}


		printModelStateMatlabFormat( q.data(), Nmu, Nnu, Nphi, NmuLid, NphiLid, lv, rm, outputName );

		double At = 0;
		printModelValuesMatlabFormat( q.data(), At, lv, rm, outputName2 );

	}

	// clean up
	lv.clearLV();

}





void printStaticLVShapeMatlabFormatFixedRot( vector <double> & q, vector <double> & theta, vector <double> & shift, string outputName, DataContainer & prmData)
{

	// initialize model
	GdmLV lv;
	lv.setupLV( prmData );


	// visualize model
	int Nmu = 10;
	int Nnu = 20;
	int Nphi = 20;
	int NmuLid = 20;
	int NphiLid = Nphi;


	RigidMotion rm;
	rm.computeRigidPrm( theta.data(), shift.data() );


	printModelStateMatlabFormat( q.data(), Nmu, Nnu, Nphi, NmuLid, NphiLid, lv, rm, outputName );

	// clean up
	lv.clearLV();

}


void printInitialLVShapeMatlabFormat( string outputName, DataContainer & prmData)
{
	// initialize model
	GdmLV lv;
	lv.setupLV( prmData );

	vector<double> q(lv.prm.Nq,0);


	// Compute the rigid motion rotation matrices
	vector <double> theta(3,0);
	vector <double> shift(3,0);

	string initSurfacePrefix = "initSurf";
	readPrmVector("theta",  initSurfacePrefix, prmData, theta);
	readPrmVector("shift",  initSurfacePrefix, prmData, shift);

	RigidMotion rm;
	rm.computeRigidPrm( theta.data(), shift.data() );


	// visualize model
	int Nmu = 10;
	int Nnu = 20;
	int Nphi = 20;
	int NmuLid = 20;
	int NphiLid = Nphi;
	printModelStateMatlabFormat( q.data(), Nmu, Nnu, Nphi, NmuLid, NphiLid, lv, rm, outputName );


	// clean up
	lv.clearLV();

}


void printModelStateMatlabFormat( double * q, int Nmu, int Nnu, int Nphi, int NmuLid, int NphiLid, GdmLV & lv, RigidMotion & rm, string outputName )
{


	if ( getProcID() == ROOT_ID )
	{
		// arrays
		double ** xEndo = new double * [Nnu];
		double ** yEndo = new double * [Nnu];
		double ** zEndo = new double * [Nnu];
		double ** xEpi = new double * [Nnu];
		double ** yEpi = new double * [Nnu];
		double ** zEpi = new double * [Nnu];
		for(int k = 0; k < Nnu; k++ )
		{
			xEndo[k] = new double [Nphi];
			yEndo[k] = new double [Nphi];
			zEndo[k] = new double [Nphi];
			xEpi[k] = new double [Nphi];
			yEpi[k] = new double [Nphi];
			zEpi[k] = new double [Nphi];
		}


		double ** xTop = new double * [Nmu];
		double ** yTop = new double * [Nmu];
		double ** zTop = new double * [Nmu];
		for(int k = 0; k < Nmu; k++ )
		{
			xTop[k] = new double [Nphi];
			yTop[k] = new double [Nphi];
			zTop[k] = new double [Nphi];
		}

		double ** xLid = new double * [NmuLid];
		double ** yLid = new double * [NmuLid];
		double ** zLid = new double * [NmuLid];
		for(int k = 0; k < NmuLid; k++ )
		{
			xLid[k] = new double [NphiLid];
			yLid[k] = new double [NphiLid];
			zLid[k] = new double [NphiLid];
		}




		FourierDeformation fourierDef;
		fourierDef.setSize( lv.prm.muinNnuGrid, lv.prm.muinNphiGrid, lv.prm.nuNnuGrid, lv.prm.nuNphiGrid, lv.prm.phiNnuGrid, lv.prm.phiNphiGrid);

		NodeGdmLV node;
		// Initialize the nodes and the static values
		node.initializeTensors( & lv.prm );



		for( int i = 0; i < Nmu; i++ )
		{
			for(int j = 0; j < Nnu; j++)
			{
				for(int k = 0; k < Nphi; k++)
				{
					node.setFixedLocation( i, j, k, Nmu, Nnu, Nphi, &lv.prm);
					node.initialComputations( &lv.prm );
					node.computeDeformation( q, fourierDef, &lv.prm );
					node.computeSubvalues();

					double x = node.x;
					double y = node.y;
					double z = node.z;

					double xRot, yRot, zRot;
					rm.computeRigid(x,y,z,xRot,yRot,zRot);


					// save endo values
					if( i == 0 )
					{
						xEndo[j][k] = xRot;
						yEndo[j][k] = yRot;
						zEndo[j][k] = zRot;
					}

					// save epi values
					if( i == Nmu-1 )
					{
						xEpi[j][k] = xRot;
						yEpi[j][k] = yRot;
						zEpi[j][k] = zRot;
					}

					// save top values
					if( j == 0 )
					{
						xTop[i][k] = xRot;
						yTop[i][k] = yRot;
						zTop[i][k] = zRot;
					}
				}
			}
		}


		// Lid export
		LidGdmLV lid;

		// create lid
		lid.setLidSize( NmuLid, NphiLid );
		lid.createLid( lv.prm );

		// set initial positions
		lid.setInitialPositions( lv.prm );


		lid.computeLid( q, fourierDef, lv.prm );

		for( int i = 0; i < lid.NmuLid; i++ )
		{
			for(int j = 0; j < lid.Nphi; j++ )
			{

				double x = lid.x[i][j];
				double y = lid.y[i][j];
				double z = lid.z[i][j];

				double xRot, yRot, zRot;
				rm.computeRigid(x,y,z,xRot,yRot,zRot);


				xLid[i][j] = xRot;
				yLid[i][j] = yRot;
				zLid[i][j] = zRot;

			}
		}



		// EXPORT

		ofstream mFile( outputName.c_str() , ios::out );


			// void printMatlab2DArray( ofstream &fileID, string key, double ** values, int N1, int N2 );
			printMatlab2DArraySimple(mFile, "xEndo", xEndo, Nnu, Nphi);
			printMatlab2DArraySimple(mFile, "yEndo", yEndo, Nnu, Nphi);
			printMatlab2DArraySimple(mFile, "zEndo", zEndo, Nnu, Nphi);

			printMatlab2DArraySimple(mFile, "xEpi", xEpi, Nnu, Nphi);
			printMatlab2DArraySimple(mFile, "yEpi", yEpi, Nnu, Nphi);
			printMatlab2DArraySimple(mFile, "zEpi", zEpi, Nnu, Nphi);

			printMatlab2DArraySimple(mFile, "xTop", xTop, Nmu, Nphi);
			printMatlab2DArraySimple(mFile, "yTop", yTop, Nmu, Nphi);
			printMatlab2DArraySimple(mFile, "zTop", zTop, Nmu, Nphi);

			printMatlab2DArraySimple(mFile, "xLid", xLid, lid.NmuLid, Nphi);
			printMatlab2DArraySimple(mFile, "yLid", yLid, lid.NmuLid, Nphi);
			printMatlab2DArraySimple(mFile, "zLid", zLid, lid.NmuLid, Nphi);

		mFile.close();








		// clean up
		node.destroyTensors();
		lid.destroyLid();

		for(int k = 0; k < Nnu; k++ )
		{
			delete [] xEndo[k];
			delete [] yEndo[k];
			delete [] zEndo[k];
			delete [] xEpi[k];
			delete [] yEpi[k];
			delete [] zEpi[k];
		}
		delete [] xEndo;
		delete [] yEndo;
		delete [] zEndo;
		delete [] xEpi;
		delete [] yEpi;
		delete [] zEpi;


		for(int k = 0; k < Nmu; k++ )
		{
			delete [] xTop[k];
			delete [] yTop[k];
			delete [] zTop[k];
		}

		delete [] xTop;
		delete [] yTop;
		delete [] zTop;




		for(int k = 0; k < Nmu; k++ )
		{
			delete [] xLid[k];
			delete [] yLid[k];
			delete [] zLid[k];
		}
		delete [] xLid;
		delete [] yLid;
		delete [] zLid;
	}

}








void printModelValuesMatlabFormat( double * q, double At, GdmLV & lv, RigidMotion & rm, string outputName )
{
	 int Nmu = lv.prm.Nmu;
	 int Nnu = lv.prm.Nnu;
	 int Nphi = lv.prm.Nphi;

	if ( getProcID() == ROOT_ID )
	{
		// arrays
		double *** x = new double ** [Nmu];
		double *** y = new double ** [Nmu];
		double *** z = new double ** [Nmu];
		double *** values1 = new double ** [Nmu];
		double *** IV0 = new double ** [Nmu];

		for(int i = 0; i < Nmu; i++ )
		{
			x[i] = new double * [Nnu];
			y[i] = new double * [Nnu];
			z[i] = new double * [Nnu];
			values1[i] = new double * [Nnu];
			IV0[i] = new double * [Nnu];

			for(int j = 0; j < Nnu; j++)
			{
				x[i][j] = new double [Nphi];
				y[i][j] = new double [Nphi];
				z[i][j] = new double [Nphi];
				values1[i][j] = new double [Nphi];
				IV0[i][j] = new double [Nphi];

			}
		}




		FourierDeformation fourierDef;
		fourierDef.setSize( lv.prm.muinNnuGrid, lv.prm.muinNphiGrid, lv.prm.nuNnuGrid, lv.prm.nuNphiGrid, lv.prm.phiNnuGrid, lv.prm.phiNphiGrid);

		NodeGdmLV node;
		// Initialize the nodes and the static values
		node.initializeTensors( & lv.prm );



		for( int i = 0; i < Nmu; i++ )
		{
			for(int j = 0; j < Nnu; j++)
			{
				for(int k = 0; k < Nphi; k++)
				{
					// set position
					node.setFixedLocation( i, j, k, Nmu, Nnu, Nphi, &lv.prm);
					node.initialComputations( &lv.prm );

					// compute values
					node.deformationStrainComputations( q, fourierDef, &lv.prm );
					node.computeStressesTractions( &lv.prm, 0, At );

					double xmodel = node.x;
					double ymodel = node.y;
					double zmodel = node.z;

					double xRot, yRot, zRot;
					rm.computeRigid(xmodel,ymodel,zmodel,xRot,yRot,zRot);


					x[i][j][k] = xRot;
					y[i][j][k] = yRot;
					z[i][j][k] = zRot;

					values1[i][j][k] = node.lambda;

					IV0[i][j][k] = node.IV0weight;

				}
			}
		}




		// EXPORT

		ofstream mFile( outputName.c_str() , ios::out );

			// void printMatlab2DArray( ofstream &fileID, string key, double ** values, int N1, int N2 );
			printMatlab3DArraySimple(mFile, "x", x, Nmu, Nnu, Nphi);
			printMatlab3DArraySimple(mFile, "y", y, Nmu, Nnu, Nphi);
			printMatlab3DArraySimple(mFile, "z", z, Nmu, Nnu, Nphi);

			printMatlab3DArraySimple(mFile, "values1",values1, Nmu, Nnu, Nphi);
			printMatlab3DArraySimple(mFile, "IV0", IV0, Nmu, Nnu, Nphi);

		mFile.close();



		// clean up
		node.destroyTensors();

		for(int i = 0; i < Nmu; i++ )
		{

			for(int j = 0; j < Nnu; j++)
			{
				delete [] x[i][j];
				delete [] y[i][j];
				delete [] z[i][j];
				delete [] values1[i][j];
				delete [] IV0[i][j];

			}

			delete [] x[i];
			delete [] y[i];
			delete [] z[i];
			delete [] values1[i];
			delete [] IV0[i];
		}
		delete [] x;
		delete [] y;
		delete [] z;
		delete [] values1;
		delete [] IV0;

	}
}
