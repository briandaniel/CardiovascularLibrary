/*
 * cdmlv_node_output.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: brian
 */

#include "cdmlv_node.hpp"





void CDMNode::printNodeToMFile( ofstream &fileID ){

	std::string structName;
	std::string key;
	structName = "nodeStruct.node" + to_string(globalID) + '.';

	// 1. indices
	key = "i";
	printMatlabVariableSimple( fileID, structName + key, imu );
	key = "j";
	printMatlabVariableSimple( fileID, structName + key, jnu );
	key = "k";
	printMatlabVariableSimple( fileID, structName + key, kphi );

	// 2. reference prolate coordinates
	key = "mu0";
	printMatlabVariableSimple( fileID, structName + key, mu0 );
	key = "nu0";
	printMatlabVariableSimple( fileID, structName + key, nu0 );
	key = "phi0";
	printMatlabVariableSimple( fileID, structName + key, phi0 );

	// 3. reference cartesian coordinates
	key = "x0";
	printMatlabVariableSimple( fileID, structName + key, x0 );
	key = "y0";
	printMatlabVariableSimple( fileID, structName + key, y0 );
	key = "z0";
	printMatlabVariableSimple( fileID, structName + key, z0 );

	// 4. deformed prolate coordinates
	key = "mu";
	printMatlabVariableSimple( fileID, structName + key, mu );
	key = "nu";
	printMatlabVariableSimple( fileID, structName + key, nu );
	key = "phi";
	printMatlabVariableSimple( fileID, structName + key, phi );

	// 5. deformed cartesian coordinates
	key = "x";
	printMatlabVariableSimple( fileID, structName + key, x );
	key = "y";
	printMatlabVariableSimple( fileID, structName + key, y );
	key = "z";
	printMatlabVariableSimple( fileID, structName + key, z );

	// 6. surface normals
	key = "nSide";
	printMatlab1DArray( fileID, structName + key, nSide );
	key = "nTop";
	printMatlab1DArray( fileID, structName + key, nTop );

	// 7. body forcing
	key = "bodyForcing";
	printMatlab1DArray( fileID, structName + key, bodyForcing);

	// 8. surface traction forcing
	key = "tractionForcingTop";
	printMatlab1DArray( fileID, structName + key, tractionForcingTop);
	key = "tractionForcingSide";
	printMatlab1DArray( fileID, structName + key, tractionForcingSide);


	// 9. cartesian representation of the divergence of the Pk1 stress
	key = "divSigma_cart";
	printMatlab1DArray( fileID, structName + key, divSigma_cart);


	// 10. displacement vector
	key = "udisp";
	printMatlab1DArray( fileID, structName + key, udisp);


	// 11. deformation gradient tensor
	key = "F";
	printMatlab2DArray( fileID, structName + key, F );


	// 12. Cauchy stress
	computeCauchyStress(); // compute cauchy stress
	key = "sigma_prol";
	printMatlab2DArray( fileID, structName + key, sigma_prol );
	key = "sigma_cart";
	printMatlab2DArray( fileID, structName + key, sigma_cart );


	// 13. fiber stretch
	key = "lambda";
	printMatlabVariableSimple( fileID, structName + key, lambda );

	// 14. green strain
	key = "Eprol";
	printMatlab2DArray( fileID, structName + key, E);

	// 15. PK2 elastic stress
	key = "Se";
	printMatlab2DArray( fileID, structName + key, Se);

}









