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
struct BurnupSim {
	double burnoutTime;//The time the fire went out = length of the fire (s?).
                       //A value of -1 indicates the fuel did not ignite.  A value
                       //of -2 indicates the fuel did not complete drying.  In such
                       //cases most of remaining return variables will not be
                       //meaningful.



	//Calculated outputs by interaction pairs:
	std::vector<double> w_o_kl;// = wo, Current (final?) ovendry loading for the larger of each component pair, kg/sq m.
	//std::vector<double> w_o_out;//The remaining ovendry loading after the fire (kg/sq m).
	std::vector<double> xmat;
	std::vector<double> tign_kl;
	std::vector<double> tout_kl;
	std::vector<double> diam;
	
	//Output by fuel type:
	std::vector<double> w_o_ij_Initial;//The amount of fuels pre-burn by fuel type (kg/m^2).
	std::vector<double> w_o_ij_Final;//The amount of fuels remaining post-burn by fuel type (kg/m^2).
	std::vector<double> combustion_ij;//The amount of fuels combusted by fuel type (kg/m^2).
	
	std::vector<double> tign_ij;//(Minimum) ignition time by fuel type (s).
	tout_ij_Min;//The minimum burnout time among all pairs within a fuel type (s).
	tout_ij_Max;//The maximum burnout time among all pairs within a fuel type (s).
	
	//Fire history to be added...
}

#endif
