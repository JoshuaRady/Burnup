/***************************************************************************************************
BurnupHistory.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/18/2025
Reference: Proj. 11 Exp. 25

	This provides an object to store Burnup model state during the time evolution of a simulation.

Licence?????
***************************************************************************************************/

#include "BurnupHistory.h"

//We need a persistant object to short data to.  This is kept private in this file, only being\
//accessed via the provided functions:
BurnupHistory BUHistStore;


//void BurnupHistory::SetTimeSteps(const int numTimeSteps)

//SaveStateToFile():
/** Store the state of the simulation at the current timestep to a BurnupHistory object for later
 * use.  Sequential calls to this routine will produce a full history of the simulated fire.  This
 * is provided as a programatic alternative to saving the history to a file with SaveStateToFile().
 *
 * @param[in] ts		Current timestep count.
 * @param[in] time		Current time (s).
 * @param[in] number	Actual number of fuel components.
 * @param[in] parts		Fuel component names / labels. [maxno]
 * @param[in] wo		Current ovendry loading for the larger of each component pair, kg / sq m. [maxkl]
 * @param[in] diam		Current diameter of the larger of each fuel component pair, m. [maxkl]	!!!!!
 * @param[in] fi		Current fire intensity (site avg), kW / sq m.
 *
 * @returns Nothing.
 */
void SaveStateToHistory(const int ts, const double time, const int number,
                        const std::vector<std::string> parts, const std::vector<double> wo,
                        const std::vector<double> diam, const double fi)
{
	/*Each call to this function stores a new time step of data to the history.  We need to make
	sure there is room for the data.  We can do this in two ways.  We could increase the size of
	the data stores each time we add a point, but that is not very efficient with vectors.
	Alternatively, we can size the data at the outset.  Since the number of timesteps is known at
	the outset of a Burnup simulations this seems like the right way to go, but it wouldn't be a bad
	idea to check as we go.*/
	
}

/** Get the history for the last simulation.
 *
 */
BurnupHistory GetHistory()
{
	//Add checking that the history is complete?
	return BUHistStore;
}
