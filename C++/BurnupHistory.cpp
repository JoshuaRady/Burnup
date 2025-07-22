/***************************************************************************************************
BurnupHistory.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/18/2025
Reference: Proj. 11 Exp. 25

	This provides an object to store Burnup model state during the time evolution of a simulation.

Licence?????
***************************************************************************************************/

#include <cmath>//For fabs().
#include "BurnupHistory.h"

//The default(ish) number of timesteps is 3000.  We add one since we currently also record the
//initial state as well.
const int NumTimeStepsDefault = 3001;

//We need a persistant object to short data to.  This is kept private in this file, only being
//accessed via the provided functions:
//The downside of the this approach is that it will sit around taking up a good bit of space, even
//when not in use.
BurnupHistory BUHistStore;

/** Default constructor
 */
BurnupHistory::BurnupHistory()
{
	//The object starts off empty but we reserve about the amount of space we expect to use:
	timestep.reserve(NumTimeStepsDefault);
	timeSec.reserve(NumTimeStepsDefault);
	fireIntensity.reserve(NumTimeStepsDefault);
}

//We could add a way to set up or reserve more space than the default.  However, tt is not clear yet
//if this is needed.
//void BurnupHistory::SetTimeSteps(const int numTimeSteps)//numFuelsTypes

/** Is the object currently empty?
 *
 * @returns Nothing.
 */
bool BurnupHistory::Empty() const
{
	return timestep.empty();
}

/** Add the state of a Burnup simulation for a timestep.  Sequential calls to this routine will
 * produce a full history of the simulated fire.
 *
 * @param[in] ts			Current timestep count.
 * @param[in] time			Current time (s).
 * @param[in] numFuelsTypes	Actual number of fuel components.		Or numFuels?????
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
                                const std::vector<std::string>& parts, const std::vector<double>& wo,
                                const double fi)
{
	/*Each call to this function stores a new timestep of data to the history.  By reserving a
	reasonable amount of space we can add length to our vectors efficiently using push_back().
	Since the number of timesteps is known at the outset of a simulation and we generally use the
	default, we do this in the constructor*/
	timestep.push_back(ts);
	timeSec.push_back(time);

	//Store the fuel loading: ToDo!!!!!
	//If this is the first timestep we will need to create vectors for each fuel type. (Or do this earlier.)
	//Possibly record the fuel names?
	//Sum the loadings for each fuel type across all components and store.

	fireIntensity.push_back(fi);
}

/** Calculate total energy produced by the fire from the fire intensity history.
 *
 * @returns The total energy released during the fire, including that of the flaming front (kJ/m^2).?????
 */
double BurnupHistory::IntegrateFireIntensity() const
{
	if (Empty())
	{
		//Report an error!!!!!
		return 0.0;
	}
	
	/*We include the igniting fire intensity from the flmaming front, which is an input to Burnup
	rather than a computed output, at the start of the history at timestep 0.  The flaming front is
	modeled as an intensity and residence time, which varies in length.  We treat the intensity for
	this timestep as constant so the cummulative energy is the intensity times the residence time,
	which we can recover from the time at timestep 1.*/
	double totalEnergy = fireIntensity[0] * timeSec[1];

	//The remaining timesteps represent Burnup calculated output and will have regular timesteps.
	//We can recover that by looking at the next two and assuming they are all the same from there:
	double dT = timeSec[2] - timeSec[1];
	
	/*The intensity tends to drop very quickly from the igniting intensity and then drop more
	slowly from there.  We could integrate this as a stepped curve, assuming the intensity is
	constant for each timestep but it is more realistic to interpolate between each point and
	integrate the area of each quadrilateral.  The only issue is what to do with the last value,
	where we have no futher value after it.  From Burnup's description the last point seems to
	represent the final intensity at burnout.  We can therefore ignore anything beyond that.*/
	for (int i = 1; i < timestep.size() - 2; i++)
	{
		double midHeight = (fireIntensity[i] + fireIntensity[i + 1]) / 2;
		totalEnergy += midHeight * dT;
	}

	return totalEnergy;
}

//External functions:-------------------------------------------------------------------------------

/** Store the state of a Burnup simulation at the current timestep to a BurnupHistory object for
 * later use.  Sequential calls to this routine will produce a full history of the simulated fire.
 * This is provided as a programatic alternative to saving the history to a file with
 * SaveStateToFile().
 * 
 * @par The level of detail stored is less than in SaveStateToFile().  Fire intensity and fuel
 * loading over time are recorded as these are the most important features needed to understand the
 * simulation.  SaveStateToFile() stores the fuel loadings and fuel diameters by component pairs at
 * each timestep. This function will only store the total loadings for each fuel class.  It doesn't
 * store diameters as there is no useful way to average diameter across pairs for a fule type.
 *
 * @param[in] ts			Current timestep count.
 * @param[in] time			Current time (s).
 * @param[in] numFuelTypes	Actual number of fuel components.
 * @param[in] parts			Fuel component names / labels. [maxno]
 * @param[in] wo			Current ovendry loading for the larger of each component pair, kg / sq m. [maxkl]
 * @param[in] fi			Current fire intensity (site avg), kW / sq m.
 *
 * @returns Nothing.
 * 
 * @note This is a public wrapper for access to the hidden private BUHistStore instantiation.
 */
void SaveStateToHistory(const int ts, const double time, const int numFuelTypes,
                        const std::vector<std::string>& parts, const std::vector<double>& wo,
                        const double fi)
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
