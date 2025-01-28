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
 * This object is used by BurnupFM() to organize and return Burnup outputs.  It also contains
 * additional processed outputs that make the data easier to use and printing utilities to display
 * the output as text. History output will be added in the future.
 *
 * @note This data structure may need some refinement.  It has been designed to deal with several
 * different use cases but it may be a bit confusing.  We will review it after using it a bit.
 */
struct BurnupSim {
	double burnoutTime;//The time the fire went out = length of the fire (s?).
	                   //A value of -1 indicates the fuel did not ignite.  A value
	                   //of -2 indicates the fuel did not complete drying.  In such
	                   //cases most of remaining return variables will not be
	                   //meaningful.
	
	double finalFireIntensity;//The final fire intensity (kW / m^2).

	//Data by fuel type:
	//The ij suffix indicate variables that are organized by fuel type.
	//This notation may be a bit confusing when the output is returned in Burnup form as the fuel
	//types may not truly be fuel model ij order.  Is there a better notation?

	bool fuelModelFormat;//Or burnupFormat?
	                     //If true the fuel level output is in the same units and fuel order as the
	                     //fuel model input to Burnup.  If false the output is in the original
	                     //Burnup units and the fuels may be reordered.

	//Select fuel properties (inputs) by type:
	//Some of the fuel properties are saved to make it easier to interpret the outputs.
	std::vector<double> SAV_ij;//Characteristic surface-area-to-volume ratios for each fuel type (cm^2/cm^3).
	std::vector<double> M_f_ij;//Fuel moisture content for each fuel type (fraction: water weight/dry fuel weight).

	std::vector<std::string> fuelNames;//Fuel type names: Burnup specific.

	//Outputs by fuel type:
	std::vector<double> w_o_ij_Initial;//The amount of fuel pre-burn by fuel type (kg/m^2).
	std::vector<double> w_o_ij_Final;//The amount of fuel remaining post-burn by fuel type (kg/m^2).
	std::vector<double> combustion_ij;//The amount of fuel combusted by fuel type (kg/m^2).
	
	std::vector<double> tign_ij;//(Minimum) ignition time by fuel type (s).
	std::vector<double> tout_ij_Min;//The minimum burnout time for all pairs within a fuel type (s).
	std::vector<double> tout_ij_Max;//The maximum burnout time for all pairs within a fuel type (s).

	//Add final diameters...

	//Calculated outputs by interaction pairs:
	//These are raw outputs from Burnup and the fuels may have been reordered.  If klFuelsReordered
	//is true use klFuelNames to interpret the fuel indexes.
	//The units for these variables are consistent between fuel models and Burnup.
	//We use the kl suffix to distinguish the fact that these variables are in kl matrix space.

	bool klFuelsReordered;//Did Burnup reorder the fuels?
	                      //This will remain true even if the fuels are reorder for the ij variables.
	                      //Do we need to specify kl?  Is there a clearer way to note this state?
	std::vector<std::string> klFuelNames;//Change to or add an order?
	//This is only needed if fuelModelFormat = true and klFuelsReordered = true.

	std::vector<double> w_o_kl;//Final ovendry loading for the larger of each component pair (kg/sq m). [wo in original Burnup]
	std::vector<double> xmat_kl;//Table of influence fractions between components.
	std::vector<double> tign_kl;//Ignition time for the larger of each fuel component pair (s).
	std::vector<double> tout_kl;//Burnout time of larger component of pairs (s).
	std::vector<double> diam_kl;//Final diameter of the larger of each fuel component pair (m).//!!!!!!

	//Fire history to be added...

	std::ostream& Print(std::ostream& output) const;
};

std::ostream& operator<<(std::ostream& output, const BurnupSim& fm);

BurnupSim BurnupFM(FuelModel fuelModel, const double duffLoading, const double duffMoisture,
                   const double tempAirC, const double U, const double fireIntensity,
                   const double t_r, const double dT, const int nTimeSteps,
                   const bool burnupFormat = false,
                   const std::vector <double> tpig_ij = {},
                   const double ak = -1.0, const double r0 = 1.83, const double dr = 0.40);

void Reorder(std::vector<double>& vec, const std::vector<int> order);

#endif
