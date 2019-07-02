/*
 * cdmlv_param.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */


#include "cdmlv_param.hpp"

void ParamCDMLV::importParams( string prmPrefix, DataContainer & prmData ){

	// Domain parameters
	Nmu = (int) readPrmValue( "Nmu", prmPrefix, prmData);
	Nnu = (int) readPrmValue( "Nnu", prmPrefix, prmData);
	Nphi = (int) readPrmValue( "Nphi", prmPrefix, prmData);
	NmuLid = (int) readPrmValue( "NmuLid", prmPrefix, prmData);

	// Pull deformation parameters
	muinNnuGrid = (int) readPrmValue( "muinNnuGrid", prmPrefix, prmData);
	muinNphiGrid = (int) readPrmValue( "muinNphiGrid", prmPrefix, prmData);
	nuNnuGrid = (int) readPrmValue( "nuNnuGrid", prmPrefix, prmData);
	nuNphiGrid = (int) readPrmValue( "nuNphiGrid", prmPrefix, prmData);
	phiNnuGrid = (int) readPrmValue( "phiNnuGrid", prmPrefix, prmData);
	phiNphiGrid = (int) readPrmValue( "phiNphiGrid", prmPrefix, prmData);
	dq = readPrmValue( "dq", prmPrefix, prmData);

	// Pull fiber params
	psi_in_b0 = readPrmValue( "psi_in_b0", prmPrefix, prmData );
	psi_out_b0 = readPrmValue( "psi_out_b0", prmPrefix, prmData );
	nuTrans = readPrmValue( "nuTrans", prmPrefix, prmData );

	// Pull stress params
	ke = readPrmValue( "ke", prmPrefix, prmData );
	bxx = readPrmValue( "bxx", prmPrefix, prmData );
	bff = readPrmValue( "bff", prmPrefix, prmData );
	bfx = readPrmValue( "bfx", prmPrefix, prmData );

	kv = readPrmValue( "kv", prmPrefix, prmData );
	ka = readPrmValue( "ka", prmPrefix, prmData );
	kav = readPrmValue( "kav", prmPrefix, prmData );

	Ls0 = readPrmValue( "Ls0", prmPrefix, prmData );
	Lsmax = readPrmValue( "Lsmax", prmPrefix, prmData );
	Lsw = readPrmValue( "Lsw", prmPrefix, prmData );



	// Read surface values from file
	surfaceInputFileName = readPrmString( "surfaceInputFileName", prmPrefix, prmData);
	DataContainer surfData;
	surfData.readDataFileMPI(surfaceInputFileName, 0);

	a = surfData.readDataValue( "a", 0);

	NtopSplines = (int) surfData.readDataValue( "NtopSplines", 0);
	Ne0 = (int) surfData.readDataValue( "Ne0", 0);
	NnuSplines = (int) surfData.readDataValue( "NnuSplines", 0);
	NphiSplines = (int) surfData.readDataValue( "NphiSplines", 0);

	muin0_const = surfData.readDataValue( "muin0_const", 0);
	muout0_const = surfData.readDataValue( "muout0_const", 0);
	muSplineConst = surfData.readDataValue( "muSplineConst", 0);
	nuUp0Min = surfData.readDataValue( "nuUp0Min", 0);

	surfData.copy1DDataSet( "e0", e0 );
	surfData.copy1DDataSet( "cin0", cin0 );
	surfData.copy1DDataSet( "cout0", cout0 );



	// These computations are required to properly setup the surfaces/deformation
	setupComputations();
}




void ParamCDMLV::setupComputations()
{
	// Set the size of the fourier deformation
	fourierDef.setSize( muinNnuGrid, muinNphiGrid, nuNnuGrid, nuNphiGrid, phiNnuGrid, phiNphiGrid );

	// Compute the number of deformation parameters
	fourierDef.computeNParams(Nq,Nc,Nd,Ne);
	if(getProcID() == ROOT_ID)
	{
		cout << "Using fourier mode deformation with Nq = " << Nq << endl;
		cout << "   with component sizes Nc, Nd, Ne = " << Nc << ", " << Nd << ", " << Ne << " modes." << endl;

	}

	// Compute top splines (this also computes nuUp0Min)
	computeTopSplineFromVector( e0, NtopSplines, nuUp0Spline, nuUp0Min );

	// generate the endo/epi surface splines
	endo0.createRegularSplines( cin0.data(), NnuSplines, NphiSplines, nuUp0Min, muSplineConst );
	epi0.createRegularSplines( cout0.data(), NnuSplines, NphiSplines, nuUp0Min, muSplineConst );

}



// Default constructor
ParamCDMLV::ParamCDMLV(){ // @suppress("Class members should be properly initialized")

	// Domain size
	Nmu = 7;
	Nnu = 11;
	Nphi = 15;
	NmuLid = 11;

	// Prolate coordinate system size
	a = 5.69194325138314;

	// Surface variables
	NnuSplines = 4;
	NphiSplines = 4;
	NtopSplines = 4;
	Ne0 = 8;

	muin0_const = 0.2998297455337583;
	muout0_const = 0.4981792755471056;
	muSplineConst = 0.2998297455337583;
	nuUp0Min = 1.056416611824186;

	// Surface definitions
	e0 = {1.203741801144519, 1.169409291041957, 1.060874629105798, 1.113358595623412, 0.06317245187642953, -0.1306440669920655, -0.03083353008982437, 0.07368285532146948};
	cin0 = {0.2678932420316608, 0.2776115477732678, -0.02837685476499962, 0.3096597167867517, 0.3700968201827635, 0.3252440788399363, 0.2650684721599101, 0.2979726177728104,
			0.3155611780275802, 0.3247720237277021, 0.3168763793958357, 0.3300702245611955, 0.3037939563649455, 0.3049101933813005, 0.3227954018475644, 0.2475198086556444,
			0.2534648859277491, 0.2664233431853664, 0.2580294444919503, 0.08211296469229715, -0.3579877520816578, 0.1031649394492068, 0.2700739674267696, 0.1266835141415793,
			-0.07422353588508093, -0.06963353996921305, 0.1518336913800593, -0.08182988121060469, -0.08143605588084647, -0.09899898054230627, -0.1618240462847307,
			-0.1369894983211486, -0.1468051828204175, -0.1743007847033025, -0.1651443479189406, 0.05290799871146508, 0.005759781466386464, -0.05319166169937021,
			-0.005100960443147073, 0.00614573800896286, 0.004536292886628992, 0.01299557540869823, -0.0251899446652754, 0.001094712620055147, -0.0255300799399596,
			0.02369453864081113, 0.0001875692650880621, 0.003948516801993364, 0.005365980307693004, 0.009458471409597227, -0.01877125539313639, -0.5428750331932296,
			0.2743854684130966, 0.1564121238776178, 0.1326778396414957, -0.07438491234836445, -0.131278431163528, 0.1416383736144125, 0.06898805167482783, 0.02258036486130348,
			0.01231327916781706, -0.0592880360323484, 0.02916175946895885, -0.0008308633903865064, -0.005953954497259682, -0.03308451137995196, 0.03817902790755524};

	cout0 = { 0.5542218334584865, 0.09474481519420352, -0.01570519946864509, 0.4531686533137138, 0.4211173259009935, 0.5130266578273666, 0.4248561016688028, 0.5204008978538351,
			0.4958110965089543, 0.4995102450697975, 0.4941733393878986, 0.5407290556618309, 0.5003849970817443, 0.4892379135847675, 0.4990557168342128, 0.4792674130733982,
			0.4826711118966646, 0.4859785219743141, 0.4608232731893663, -0.02135756940425587, 0.1028789939983993, 0.06326933482311128, 0.1908415624130073, 0.1867904019975304,
			0.100296912451799, -0.1743712565775788, 0.1663974226832234, -0.09114215235716895, -0.03172029914598919, 0.05988233981737087, -0.1840434551502231, -0.02063871006112681,
			0.0405192908929534, 0.08323881677651829, 0.03231126709654744, -0.121220286903152, 0.04793523094345736, 0.03830933154555893, -0.0470959183212255, -0.0418783452457762,
			-0.04995770819304787, 0.05689279534402691, 0.006535721443067747, -0.04818932218299128, -0.05357974229060824, 0.0358065315035791, 0.03149754514662095, -0.04643593208895178,
			0.00724725030052624, -0.05913456949923727, -0.008004959325135366, -0.0307015830440325, 0.09599049029227327, 0.005981227516188659, -0.03742244498694159, -0.03725172831409849,
			-0.2040102683517043, 0.06448045138855768, 0.1271274072337597, 0.024322992118683, 0.151248537637118, -0.1241998635590508, -0.0326523417163773, 0.003299634280903639,
			0.0398128122092458, -0.04365615735890886, -0.05645171758641145 };


	// Fiber direction parameters
	nuTrans = 1.8416;
	psi_in_b0 = 1.0472;
	psi_out_b0 = -1.0472;

	// Fourier deformation size
	muinNnuGrid = 1;
	muinNphiGrid = 1;
	nuNnuGrid = 1;
	nuNphiGrid = 1;
	phiNnuGrid = 1;
	phiNphiGrid = 0;
	dq = 1e-6;

	// Stress parameters
	ke = 1.75;
	bff = 3;
	bxx = 2;
	bfx = 1;
	kv = 0.1;
	ka = 130;
	kav = 10;

	// Muscle fiber sarcomere parameters
	Ls0 = 1.82;
	Lsmax = 2.4;
	Lsw = 0.435;


	// These computations are required to properly setup the surfaces/deformation
	setupComputations();

}







void computeTopSplineFromVector( vector<double> e0, int NtopSplines, Spline1D & nuUp0Spline, double & nuUp0Min )
{
	// this is really phi
	double * x = new double [NtopSplines + 1];
	for(int k = 0; k < NtopSplines+1; k++){ x[k] = 2*PI * (double) k /  (double) NtopSplines; }


	double * y = new double [NtopSplines + 1];
	double * dydx = new double [NtopSplines + 1];


	// first set of values are function values
	int i = 0;
	for(int k = 0; k < NtopSplines; k++){ y[k] = e0[i]; i++; }

	// second set of values are derivatives
	for(int k = 0; k < NtopSplines; k++){ dydx[k] = e0[i]; i++; }

	// periodicity
	y[NtopSplines] = y[0];
	dydx[NtopSplines] = dydx[0];

	// now compute the spline coefficients
	nuUp0Spline.computeSplineCoef( NtopSplines + 1, x, y, dydx );


	// nuUp0Min is computed by finding the smallest value of the spline
	double xMin, xMax, pMin, pMax;
	nuUp0Spline.findSplineMinimumAndMaximum(xMin,xMax,pMin,pMax);
	nuUp0Min = pMin;


	// clean up
	delete [] x;
	delete [] y;
	delete [] dydx;


}


