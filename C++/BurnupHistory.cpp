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

const int NumTimeStepsDefault = 3000;

//We need a persistant object to short data to.  This is kept private in this file, only being\
//accessed via the provided functions:
BurnupHistory BUHistStore;

/** Default constructor
 */
BurnupHistory::BurnupHistory()
{
	//The object starts off empty but we reserve about the amount of space:
	timestep.reserve(NumTimeStepsDefault);
	timeSec.reserve(NumTimeStepsDefault);
	fireIntensity.reserve(NumTimeStepsDefault);
}

//We could add a way to set up or reserve more space than the default.  However, tt is not clear yet
//if this is needed.
//void BurnupHistory::SetTimeSteps(const int numTimeSteps)

/** Add the state of a Burnup simulation at timestep.  Sequential calls to this routine will
 * produce a full history of the simulated fire.
 *
 * @param[in] ts			Current timestep count.
 * @param[in] time			Current time (s).
 * @param[in] numFuelsTypes	Actual number of fuel components.		Or numFuels
 * @param[in] parts			Fuel component names / labels. [maxno]
 * @param[in] wo			Current ovendry loading for the larger of each component pair, kg / sq m. [maxkl]
 * @param[in] fi			Current fire intensity (site avg), kW / sq m.
 *
 * @returns Nothing.
 * 
 * @note The function is not yet complete.  The fuel loading in not actually stored yet and will be
 * added in future.
 */
void BurnupHistory::AddTimeStep(const int ts, const double time, const int numFuelsTypes,
                                const std::vector<std::string> parts, const std::vector<double> wo,
                                const double fi)
{
	/*Each call to this function stores a new time step of data to the history.  We need to make
	sure there is room for the data.  We can do this in two ways.  We could increase the size of
	the data stores each time we add a point, but that is not very efficient with vectors.
	Alternatively, we can size the data at the outset.  Since the number of timesteps is known at
	the outset of a Burnup simulations this seems like the right way to go, but it wouldn't be a bad
	idea to check as we go.
	Actually, by reserving a reasonable amount of space we can add length to our vectors efficiently
	and use push_back()*/

	
}

//External functions:-------------------------------------------------------------------------------

//SaveStateToFile():
/** Store the state of a Burnup simulation at the current timestep to a BurnupHistory object for
 * later use.  Sequential calls to this routine will produce a full history of the simulated fire.
 * This is provided as a programatic alternative to saving the history to a file with
 * SaveStateToFile().
 * 
 * @par The level of detail stored is less than in SaveStateToFile().  Fire intensity and fuel
 * loading over time are recorded as these are the most important features needed to understand the
 * simulation.  SaveStateToFile() stores the fuel loadings and fuel diameters by component pairs at
 * each time step. This function will only store the total loadings for each fuel class.  It doesn't
 * store diameters as there is no useful way to average diameter across pairs for a fule type.
 *
 * @param[in] ts			Current timestep count.
 * @param[in] time			Current time (s).
 * @param[in] numFuelsTypes	Actual number of fuel components.
 * @param[in] parts			Fuel component names / labels. [maxno]
 * @param[in] wo			Current ovendry loading for the larger of each component pair, kg / sq m. [maxkl]
 * @param[in] fi			Current fire intensity (site avg), kW / sq m.
 *
 * @returns Nothing.
 * 
 * @note This is a public wrapper for access to the hidden private BUHistStore instantiation.
 */
void SaveStateToHistory(const int ts, const double time, const int numFuelsTypes,
                        const std::vector<std::string> parts, const std::vector<double> wo,
                        const std::vector<double> diam, const double fi)
{
	BUHistStore.AddTimeStep(ts, time, numFuelTypes, parts, wo, fi);
}

/** Get the history for the last simulation.
 *
 */
BurnupHistory GetHistory()
{
	//Add checking that the history is complete?
	return BUHistStore;
}
