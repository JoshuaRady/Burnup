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
 *
 * @brief A data structure that holds the output of a Burnup simulation (and a selection of inputs).
 *
 */
struct BurnupSim {
	double burnoutTime;//The time the fire went out = length of the fire (s?).
	                   //A value of -1 indicates the fuel did not ignite.  A value
	                   //of -2 indicates the fuel did not complete drying.  In such
	                   //cases most of remaining return variables will not be
	                   //meaningful.
	
	double finalFireIntensity;//The final fire intensity (kW / m^2).

	//Calculated outputs by interaction pairs:
	//We use the kl suffix to distinguish the fact that these variables are in kl matrix space.
	std::vector<double> w_o_kl;//Final ovendry loading for the larger of each component pair (kg/sq m). [wo in original Burnup]
	std::vector<double> xmat_kl;//Table of influence fractions between components.
	std::vector<double> tign_kl;//Ignition time for the larger of each fuel component pair.
	std::vector<double> tout_kl;//Burnout time of larger component of pairs.
	std::vector<double> diam_kl;//Final diameter of the larger of each fuel component pair (m).//!!!!!!

	//Output by fuel type:
	//The ij suffix indicates these variables are in same vector representation used in fuel models.

	//Fuel type names...

	std::vector<double> w_o_ij_Initial;//The amount of fuel pre-burn by fuel type (kg/m^2).
	std::vector<double> w_o_ij_Final;//The amount of fuel remaining post-burn by fuel type (kg/m^2).
	std::vector<double> combustion_ij;//The amount of fuel combusted by fuel type (kg/m^2).
	
	std::vector<double> tign_ij;//(Minimum) ignition time by fuel type (s).
	std::vector<double> tout_ij_Min;//The minimum burnout time for all pairs within a fuel type (s).
	std::vector<double> tout_ij_Max;//The maximum burnout time for all pairs within a fuel type (s).

	//Final diameters...
	
	//Fire history to be added...

	std::ostream& Print(std::ostream& output) const;
};

std::ostream& operator<<(std::ostream& output, const BurnupSim& fm);

BurnupSim SimulateFM(FuelModel fuelModel, const double duffLoading, const double duffMoisture,
                     const double tempAirC, const double U, const double fireIntensity,
                     const double t_r, const double dT, const int nTimeSteps,
                     const double ak = -1.0, const double r0 = 1.83, const double dr = 0.40);

#endif
