/*
 * gdmlv_params.cpp
 *
 *  Created on: Sep 11, 2018
 *      Author: brian
 */


#include "gdmlv_params.hpp"



void ParamGdmLV::importParams( DataContainer & prmData ){

	string prmPrefix;

	// 1----------------> First parameter file prefix
	prmPrefix = "gdmlv";

	deformationModelIndex = (int) readPrmValue( "deformationModelIndex", prmPrefix, prmData);
	LVSolutionSetting = (int) readPrmValue( "LVSolutionSetting", prmPrefix, prmData);
	updateInitialValues = (int) readPrmValue( "updateInitialValues", prmPrefix, prmData);
	plotOutput = (int) readPrmValue( "plotOutput", prmPrefix, prmData);
	useRVPressure = (int) readPrmValue( "useRVPressure", prmPrefix, prmData);

	// determines if the reference shape adjustment should be used
	useRefAdj = readPrmValue( "useReferenceShapeAdjustment", prmPrefix, prmData );
	refAdjust_muinNnuBasis = readPrmValue( "refAdjust_muinNnuBasis", prmPrefix, prmData );
	refAdjust_muinNphiBasis = readPrmValue( "refAdjust_muinNphiBasis", prmPrefix, prmData );

	vector<double> q0s; // initially the reference adjustment is set to zero
	readPrmVector(  "q0s", prmPrefix, prmData, q0s );


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

	// dummy used to automatically compute Nc, Nd, Ne
	FourierDeformation fourierTemp( muinNnuGrid, muinNphiGrid, nuNnuGrid, nuNphiGrid,
			phiNnuGrid, phiNphiGrid );


	// 666 666 666 temporary
	if ( deformationModelIndex == 1 )
	{

		fourierTemp.computeNParams(Nq,Nc,Nd,Ne);

		if(getProcID() == ROOT_ID)
		{
			cout << "Using fourier mode deformation with Nc, Nd, Ne = " << Nc << ", " << Nd << ", " << Ne << " modes." << endl;
			fourierTemp.printNModesToConsole();
		}
	}
	else
	{
		Nc = 4;
		Nd = 1;
		Ne = 2;

		// Total number of parameters
		Nq = Nc + Nd + Ne;

		if(getProcID() == ROOT_ID)
		{
			cout << "Using 7 mode deformation with Nc, Nd, Ne = " << Nc << ", " << Nd << ", " << Ne << " modes." << endl;
		}
	}





	newtMaxIter = (int) readPrmValue( "newtMaxIter", prmPrefix, prmData );


	// doubles
	Tc = readPrmValue( "Tc", prmPrefix, prmData );
	Ta = readPrmValue( "Ta", prmPrefix, prmData );

	Ncyc = readPrmValue( "Ncyc", prmPrefix, prmData );

	dtMin = readPrmValue( "dtMin", prmPrefix, prmData );
	dtMax = readPrmValue( "dtMax", prmPrefix, prmData );


	// Pull fiber params
	psi_in_b0 = readPrmValue( "psi_in_b0", prmPrefix, prmData );
	psi_out_b0 = readPrmValue( "psi_out_b0", prmPrefix, prmData );
	nuTrans = readPrmValue( "nuTrans", prmPrefix, prmData );

	// Pull stress params
	ke = readPrmValue( "ke", prmPrefix, prmData );
	bxx = readPrmValue( "bxx", prmPrefix, prmData );
	bff = readPrmValue( "bff", prmPrefix, prmData );
	bfx = readPrmValue( "bfx", prmPrefix, prmData );
	rivlinMooneyMaterial = (int) readPrmValue( "rivlinMooneyMaterial", prmPrefix, prmData );

	kv = readPrmValue( "kv", prmPrefix, prmData );

	ka = readPrmValue( "ka", prmPrefix, prmData );
	kav = readPrmValue( "kav", prmPrefix, prmData );
	ked = readPrmValue( "ked", prmPrefix, prmData );
	kd = readPrmValue( "kd", prmPrefix, prmData );

	activePower = readPrmValue( "activePower", prmPrefix, prmData );


	Ls0 = readPrmValue( "Ls0", prmPrefix, prmData );
	Lsmax = readPrmValue( "Lsmax", prmPrefix, prmData );
	Lsw = readPrmValue( "Lsw", prmPrefix, prmData );

	Rmvo = readPrmValue( "Rmvo", prmPrefix, prmData );
	Rmvc = readPrmValue( "Rmvc", prmPrefix, prmData );
	Raovo = readPrmValue( "Raovo", prmPrefix, prmData );
	Raovc = readPrmValue( "Raovc", prmPrefix, prmData );
	beta = readPrmValue( "beta", prmPrefix, prmData );

	Rpao = readPrmValue( "Rpao", prmPrefix, prmData );

	deltaq = readPrmValue( "deltaq", prmPrefix, prmData );
	newtXConvg = readPrmValue( "newtXConvg", prmPrefix, prmData );
	newtHConvg = readPrmValue( "newtHConvg", prmPrefix, prmData );


	// only need output strings on root id
	matOut = readPrmString( "matOut", prmPrefix, prmData);
	plotCyclesFile = readPrmString( "plotCyclesFile", prmPrefix, prmData);

	initialValuesFile = readPrmString( "initialValuesFile", prmPrefix, prmData);
	cartesianOutFile = readPrmString( "cartesianOutFile", prmPrefix, prmData);
	matlabParamFile = readPrmString( "matlabParamFile", prmPrefix, prmData);
	staticFile = readPrmString( "staticFile", prmPrefix, prmData);
	optimParamsFile = readPrmString( "optimParamsFile", prmPrefix, prmData);
	optimOutputPrefix = readPrmString( "optimOutputPrefix", prmPrefix, prmData);
	optimOutputPostfix = readPrmString( "optimOutputPostfix", prmPrefix, prmData);

	// 2.----------------> second parameter file prefix
	prmPrefix = "lump";


	// number of lump models
	Nw = (int) readPrmValue( "Nw", prmPrefix, prmData );

	// fixed pressures (if used)
	Psa_fixed = readPrmValue( "Psa_fixed", prmPrefix, prmData );
	Ppv_fixed= readPrmValue( "Ppv_fixed", prmPrefix, prmData );

	// capacitance
	Cla = readPrmValue( "Cla", prmPrefix, prmData );
	Cpv = readPrmValue( "Cpv", prmPrefix, prmData );
	Cra = readPrmValue( "Cra", prmPrefix, prmData );
	Csa = readPrmValue( "Csa", prmPrefix, prmData );
	Csp = readPrmValue( "Csp", prmPrefix, prmData );
	Csv = readPrmValue( "Csv", prmPrefix, prmData );
	Cpp = readPrmValue( "Cpp", prmPrefix, prmData );
	Cpa = readPrmValue( "Cpa", prmPrefix, prmData );

	// unstressed volumes
	Vu_sa = readPrmValue( "Vu_sa", prmPrefix, prmData );
	Vu_sp = readPrmValue( "Vu_sp", prmPrefix, prmData );
	Vu_sv = readPrmValue( "Vu_sv", prmPrefix, prmData );
	Vu_pa = readPrmValue( "Vu_pa", prmPrefix, prmData );
	Vu_pp = readPrmValue( "Vu_pp", prmPrefix, prmData );
	Vu_pv = readPrmValue( "Vu_pv", prmPrefix, prmData );
	Vu_ra = readPrmValue( "Vu_ra", prmPrefix, prmData );
	Vu_rv = readPrmValue( "Vu_rv", prmPrefix, prmData );
	Vu_la = readPrmValue( "Vu_la", prmPrefix, prmData );

	// resistance
	Rsa = readPrmValue( "Rsa", prmPrefix, prmData );
	Rsp = readPrmValue( "Rsp", prmPrefix, prmData );
	Rsv = readPrmValue( "Rsv", prmPrefix, prmData );
	Rpa = readPrmValue( "Rpa", prmPrefix, prmData );
	Rpp = readPrmValue( "Rpp", prmPrefix, prmData );
	Rpv = readPrmValue( "Rpv", prmPrefix, prmData );
	Rra = readPrmValue( "Rra", prmPrefix, prmData );

	// right ventricle
	P0_rv = readPrmValue( "P0_rv", prmPrefix, prmData );
	ke_rv = readPrmValue( "ke_rv", prmPrefix, prmData );
	Emax_rv = readPrmValue( "Emax_rv", prmPrefix, prmData );
	kr_rv = readPrmValue( "kr_rv", prmPrefix, prmData );


	// total volume
	Vt0 = readPrmValue( "Vt0", prmPrefix, prmData );



	// 3.----------------> third parameter file prefix
	prmPrefix = "initSurf";


	NtopSplines = (int) readPrmValue( "NtopSplines", prmPrefix, prmData );
	Ne0 = 2*NtopSplines; // twice as many parameters as splines in a 1D periodic spline

	NnuSplines = (int) readPrmValue( "NnuSplines", prmPrefix, prmData );
	NphiSplines = (int) readPrmValue( "NphiSplines", prmPrefix, prmData );

	muin0_const = readPrmValue( "muin0_const", prmPrefix, prmData );
	muout0_const = readPrmValue( "muout0_const", prmPrefix, prmData );

	muSplineConst = readPrmValue( "muSplineConst", prmPrefix, prmData );

	// load value of a

	if(getProcID() == ROOT_ID)
		cout << "Finding value of a, either named a or a0." << endl;

	a = readPrmValue( "a", prmPrefix, prmData );
	double a0 = readPrmValue( "a0", prmPrefix, prmData );
	if ( a == 0 ){ a = a0;  }
	a_original = a;

	if(getProcID() == ROOT_ID)
		cout << "Set a = a0 = " << a << endl;

	theta0 = new double [3];
	shift0 = new double [3];

	readPrmVector( "theta", prmPrefix, prmData, theta0, 3 );
	readPrmVector( "shift", prmPrefix, prmData, shift0, 3 );


	// 4.----------------> fourth parameter file prefix
	prmPrefix = "staticLV";

	At_static = readPrmValue( "At_static", prmPrefix, prmData );
	Plv_static = readPrmValue( "Plv_static", prmPrefix, prmData );
	Prv_static = readPrmValue( "Prv_static", prmPrefix, prmData );


	staticOutputFileName = readPrmString( "staticOutputFileName", prmPrefix, prmData);





	// Initial computations

	twoLswSq = 2*pow(Lsw,2);


	// Spline pre-computations
	NSplines = NnuSplines*NphiSplines;
	NSplineParams = 3 + 4*NnuSplines*NphiSplines;
	Ncin0 = NSplineParams;
	Ncout0 = NSplineParams;

	e0 = new double[Ne0];
	cin0 = new double[Ncin0];
	cout0 = new double[Ncout0];

	w0 = new double [Nw];
	q0 = new double [Nq];
	for(int k = 0; k < Nq; k++){ q0[k] = 0.0;}
	dwdt0 = new double [Nw];
	dqdt0 = new double [Nq];

	NrvBndCoef = 4;
	rvBoundaryCoef = new double [NrvBndCoef];

	if(getProcID() == ROOT_ID)
	{
		cout << endl << "Number of modes ..." << endl;
		cout << "Nc = " << Nc << ", Nd = " << Nd << ",  Ne = " << Ne << ", Nq = " << Nq << endl;
		cout << "Initial values..." << endl;
		print1DArrayLine(q0,Nq,5,"q0");
	}

	// load arrays

	// 1----------------> First parameter file prefix (used again)
	prmPrefix = "gdmlv";
	readPrmVector("rvBoundaryCoef", prmPrefix, prmData, rvBoundaryCoef, NrvBndCoef );


	// 3.----------------> third parameter file prefix (used again)
	prmPrefix = "initSurf";
	readPrmVector("e0", prmPrefix, prmData, e0,  Ne0 );
	readPrmVector("cin0", prmPrefix, prmData, cin0, Ncin0 );
	readPrmVector("cout0", prmPrefix, prmData, cout0, Ncout0 );


	// 5. ----------------> fifth parameter file prefix
	prmPrefix = "initValsGDMLV";

	readPrmVector("w0", prmPrefix, prmData, w0, Nw );
	readPrmVector("q0", prmPrefix, prmData, q0, Nq );
	readPrmVector("dwdt0", prmPrefix, prmData, dwdt0, Nw );
	readPrmVector("dqdt0", prmPrefix, prmData, dqdt0, Nq );

	Plv0 = readPrmValue( "Plv0", prmPrefix, prmData );
	Ppao0 = readPrmValue( "Plv0", prmPrefix, prmData );


	// additional computations

	N = Nmu*Nnu*Nphi;


	// set top spline
	computeTopSpline();

	// setup the endo/epi surfaces

	vector<double> cin0Vec( Ncin0, 0);
	vector<double> cout0Vec ( Ncout0, 0);
	for(int i = 0; i < Ncin0; i++) { cin0Vec[i] = cin0[i]; }
	for(int i = 0; i < Ncout0; i++) { cout0Vec[i] = cout0[i]; }

	lvRef.setupRegularRefSurfaces( cin0Vec, cout0Vec, NnuSplines, NphiSplines, nuUp0Min, muSplineConst, a_original );

	// If the reference adjustment should be used, then set it up
	if( useRefAdj == 1)
	{
		lvRef.setupRefSurfAdj( refAdjust_muinNnuBasis, refAdjust_muinNphiBasis );
		int Nq0sModes = lvRef.getRefSurfAdjSize();
		q0s.resize(Nq0sModes);
		lvRef.setRefSurfAdjValues(q0s);
		lvRef.refAdjOn();
		a = lvRef.focalLength(); // set the new focal length

		if(getProcID() == ROOT_ID)
		{
			cout << "ref shape parameter is "; print1DVector(lvRef.q0s);
			cout << "Adjusted a to " << a << " from a_original = " << a_original << endl;
		}

	}
	else
	{
		lvRef.refAdjOff();
	}


	// mue is set to the size of the axisymmetric outer surface
	mue = muout0_const;



	// endo0.createConstantProlateSplines( muin0_const, NnuSplines, NphiSplines, nuUp0Min, muSplineConst );
	// epi0.createConstantProlateSplines( muout0_const, NnuSplines, NphiSplines, nuUp0Min, muSplineConst );




	// printParams();
}




void ParamGdmLV::computeTopSpline()
{

	computeTopSplineFromVector( e0, NtopSplines, &nuUp0Spline, nuUp0Min );

}


void ParamGdmLV::printMatlabParamFile(){

	ofstream outFileID( matlabParamFile.c_str() , ios::out );


	printMatlabVariable( outFileID, "useRVPressure", useRVPressure);
	printMatlabVariable( outFileID, "Nmu", Nmu );
	printMatlabVariable( outFileID, "Nnu", Nnu );
	printMatlabVariable( outFileID, "Nphi", Nphi );
	printMatlabVariable( outFileID, "NmuLid", NmuLid );
	printMatlabVariable( outFileID, "Nc", Nc );
	printMatlabVariable( outFileID, "Nd", Nd );
	printMatlabVariable( outFileID, "Ne", Ne );
	printMatlabVariable( outFileID, "Nq", Nq );

	printMatlabVariable( outFileID, "newtMaxIter", newtMaxIter);
	printMatlabVariable( outFileID, "Tc", Tc );
	printMatlabVariable( outFileID, "Ta", Ta );
	printMatlabVariable( outFileID, "Ncyc", Ncyc);
	printMatlabVariable( outFileID, "dtMin", dtMin);
	printMatlabVariable( outFileID, "dtMax", dtMax);
	printMatlabVariable( outFileID, "psi_in_b0", psi_in_b0 );
	printMatlabVariable( outFileID, "psi_out_b0", psi_out_b0 );
	printMatlabVariable( outFileID, "nuTrans", nuTrans );
	printMatlabVariable( outFileID, "ke", ke );
	printMatlabVariable( outFileID, "bxx", bxx );
	printMatlabVariable( outFileID, "bff", bff );
	printMatlabVariable( outFileID, "bfx", bfx );
	printMatlabVariable( outFileID, "kv", kv );
	printMatlabVariable( outFileID, "ka", ka );
	printMatlabVariable( outFileID, "kav", kav );
	printMatlabVariable( outFileID, "ked", ked );
	printMatlabVariable( outFileID, "kd", kd );
	printMatlabVariable( outFileID, "activePower", activePower );
	printMatlabVariable( outFileID, "Ls0", Ls0 );
	printMatlabVariable( outFileID, "Lsmax", Lsmax);
	printMatlabVariable( outFileID, "Lsw", Lsw );
	printMatlabVariable( outFileID, "Rmvo", Rmvo );
	printMatlabVariable( outFileID, "Rmvc", Rmvc);
	printMatlabVariable( outFileID, "Raovo", Raovo);
	printMatlabVariable( outFileID, "Raovc", Raovc);
	printMatlabVariable( outFileID, "beta", beta );
	printMatlabVariable( outFileID, "Rpao", Rpao );
	printMatlabVariable( outFileID, "Nw", Nw);

	printMatlabVariable( outFileID, "Psa_fixed", Psa_fixed );
	printMatlabVariable( outFileID, "Ppv_fixed", Ppv_fixed );

	printMatlabVariable( outFileID, "Cla", Cla );
	printMatlabVariable( outFileID, "Cpv", Cpv );
	printMatlabVariable( outFileID, "Cra", Cra );
	printMatlabVariable( outFileID, "Csa", Csa );
	printMatlabVariable( outFileID, "Csp", Csp );
	printMatlabVariable( outFileID, "Csv", Csv );
	printMatlabVariable( outFileID, "Cpp", Cpp );
	printMatlabVariable( outFileID, "Cpa", Cpa );

	printMatlabVariable( outFileID, "Vu_sa", Vu_sa );
	printMatlabVariable( outFileID, "Vu_sp", Vu_sp );
	printMatlabVariable( outFileID, "Vu_sv", Vu_sv );
	printMatlabVariable( outFileID, "Vu_pa", Vu_pa );
	printMatlabVariable( outFileID, "Vu_pp", Vu_pp );
	printMatlabVariable( outFileID, "Vu_pv", Vu_pv );
	printMatlabVariable( outFileID, "Vu_ra", Vu_ra );
	printMatlabVariable( outFileID, "Vu_rv", Vu_rv );
	printMatlabVariable( outFileID, "Vu_la", Vu_la);

	printMatlabVariable( outFileID, "Rsa", Rsa );
	printMatlabVariable( outFileID, "Rsp", Rsp );
	printMatlabVariable( outFileID, "Rsv", Rsv );
	printMatlabVariable( outFileID, "Rpa", Rpa );
	printMatlabVariable( outFileID, "Rpp", Rpp );
	printMatlabVariable( outFileID, "Rpv", Rpv );
	printMatlabVariable( outFileID, "Rra", Rra );



	printMatlabVariable( outFileID, "Vt0", Vt0 );
	printMatlabVariable( outFileID, "P0_rv", P0_rv );
	printMatlabVariable( outFileID, "ke_rv", ke_rv );
	printMatlabVariable( outFileID, "Emax_rv", Emax_rv );
	printMatlabVariable( outFileID, "kr_rv", kr_rv );
	printMatlabVariable( outFileID, "Ne0", Ne0 );
	printMatlabVariable( outFileID, "NnuSplines", NnuSplines );
	printMatlabVariable( outFileID, "NphiSplines", NphiSplines );
	printMatlabVariable( outFileID, "muin0_const", muin0_const );
	printMatlabVariable( outFileID, "muout0_const", muout0_const );
	printMatlabVariable( outFileID, "muSplineConst", muSplineConst );
	printMatlabVariable( outFileID, "a", a );

	printMatlabArraySimple( outFileID, "e0", e0, Ne0 );

	outFileID.close();



}



void computeTopSplineFromVector( double * e0, int NtopSplines, Spline1D * nuUp0Spline, double & nuUp0Min )
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
	nuUp0Spline->computeSplineCoef( NtopSplines + 1, x, y, dydx );


	// nuUp0Min is computed by finding the smallest value of the spline
	double xMin, xMax, pMin, pMax;
	nuUp0Spline->findSplineMinimumAndMaximum(xMin,xMax,pMin,pMax);
	nuUp0Min = pMin;


	// clean up
	delete [] x;
	delete [] y;
	delete [] dydx;


}





