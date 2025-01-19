/***************************************************************************************************
BurnupFuelModelInterface.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 1/17/2025
Reference: Proj. 11 Exp. 22

	This provides an alternate interface to the Burnup wildfire fuel consumption model that uses
fire behavior fuel models, as implemented in the Fireweed Wildfire Code Library.

Licence?????
***************************************************************************************************/
#ifndef BURNUPFUELMODELINTERFACE_H
#define BURNUPFUELMODELINTERFACE_H

#include <string>

#include "FireweedFuelModels.h"

//This will need to be moved to the header when done:
/** @struct BUSim
 * @brief A data structure that holds the output of a Burnup simulation (and a selection of inputs maybe?).
 *
 */
struct BUSim {

	//Calculated outputs
	std::vector<double> w_o_kl;// = wo, Current (final?) ovendry loading for the larger of each component pair, kg/sq m.
	//std::vector<double> w_o_out;//The remaining ovendry loading after the fire (kg/sq m).
	std::vector<double> xmat;
	std::vector<double> tign;
	std::vector<double> tout;
	std::vector<double> diam;
}


#endif
