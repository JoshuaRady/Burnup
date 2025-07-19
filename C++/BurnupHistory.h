/***************************************************************************************************
BurnupHistory.h
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/18/2025
Reference: Proj. 11 Exp. 25

	This provides an object to store Burnup model state during the time evolution of a simulation.

Licence?????
***************************************************************************************************/
#ifndef BURNUPHISTORY_H
#define BURNUPHISTORY_H

#include <vector>

/** @struct BurnupHistory
 *
 * @brief A data structure that can be used to accumulate and store the state of the Burnup model
 * simulation over time producing a full history of the simulated fire.
 *
 * This is provided as a programatic alternative to saving the history to a file with
 * SaveStateToFile().
 */
struct BurnupHistory {
	std::vector<int> timestep;//Needed?
	std::vector<double> timeSec;//The time point in the simulation (seconds).
	std::vector<double> fireIntensity;//The fire intensity (site avg) at each time step (kW/m^2).
	//parts
	//std::vector<double> w_o_ij://Vector of vectors?
	
	BurnupHistory();
	//void SetTimeSteps(const int numTimeSteps);
	void AddTimeStep(const int ts, const double time, const int numFuelTypes,
	                 const std::vector<std::string> parts, const std::vector<double> wo,
	                 const double fi);
};

//External functions:
void SaveStateToHistory(const int ts, const double time, const int numFuelTypes,
                        const std::vector<std::string> parts, const std::vector<double> wo,
                        const double fi);
BurnupHistory GetHistory();

#endif //BURNUPHISTORY_H
