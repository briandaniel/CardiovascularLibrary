/*
 * coupledSingleFiberAortaSystem.hpp
 *
 *  Created on: Jan 8, 2019
 *      Author: brian
 */

#ifndef COUPLEDCARDIACVESSEL1DSYSTEMS_SINGLEFIBERANDAORTASYSTEM_COUPLEDSINGLEFIBERAORTASYSTEM_HPP_
#define COUPLEDCARDIACVESSEL1DSYSTEMS_SINGLEFIBERANDAORTASYSTEM_COUPLEDSINGLEFIBERAORTASYSTEM_HPP_



// Parallel commands
#define ROOT_ID 0	 // id of root processor

// other includes
#include <mpi.h>
#include <cmath>

// nonlinear solver include
#include "NonlinearSolvers/NewtonsMethod.hpp"
#include "NonlinearSolvers/Nonl_Objective.hpp"

// local includes
#include "../../ActivationFunctions/temporalActivationFunctions.hpp"
#include "../../LumpModels/LumpModel.hpp"
#include "../../SingleFiberModel/SingleFiberModel.hpp"
#include "../../ValveFunctions/ValveFunctions.hpp"
#include "../../WavePropagationModel1D/WavePropagationFunctions.hpp"


// Utility library includes
#include "UtilityFunctions/utilityFunctions.hpp"
#include "ImportExport/importExport.hpp"
#include "ParameterParser/ParameterParser.hpp"

// local includes
void computeCoupledSingleFiberAortaSystem( string prmPrefix, DataContainer & prmData );


#endif /* COUPLEDCARDIACVESSEL1DSYSTEMS_SINGLEFIBERANDAORTASYSTEM_COUPLEDSINGLEFIBERAORTASYSTEM_HPP_ */
