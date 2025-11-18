/***************************************************************************************************
BurnupHistory.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/18/2025
Reference: Proj. 11 Exp. 25

	This provides an object to store Burnup model state during the time evolution of a simulation
to create a programmatically available simulation history.

Licence?????
***************************************************************************************************/

#include "BurnupHistory.h"
#include "FireweedMessaging.h"

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
 * @returns The total energy released during the fire, including that of the flaming front (kJ/m^2).
 *
 * @note The fact that we include the energy of the flaming front may be double counting.  If we
 * consider that the energy input from the flames next to the site we should also consider that a
 * similar amount of energy is lost to the adjacent patch.
 */
double BurnupHistory::IntegrateFireIntensity() const
{
	if (Empty())
	{
		//If history has not been stored return 0 and a warning.  This might warrant an error:
		Warning("BurnupHistory: The fire history is not stored.");
		return 0.0;
	}
	else if (timestep.size() == 1)
	{
		//If there is only one timestep (representing the flaming fron) the fire did not ignite:
		//We don't know why the fire didn't start without access to the parent object's burnoutTime
		//member.  The fact tha
		Warning("BurnupHistory: The fire did not ignite.");
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

/** Print the fire history data to an output stream formatted for screen reading.
 *
 * @param[in] output The output stream to print to.
 *
 * @returns The ostream so it can be concatenated to.
 */
std::ostream& BurnupHistory::Print(std::ostream& output) const
{
	output << "Burnup fire intensity history:" << std::endl;

	//Print layer properties in table form:
	const int timestepWidth = 9;//Name & description
	const int timeSecWidth = 8;//Name & units
	const int fireIntensityWidth = 14;//Name

	//Member name header:
	output << std::setw(timestepWidth) << "timestep"
		<< std::setw(timeSecWidth) << "timeSec"
		<< std::setw(fireIntensityWidth) << "fireIntensity" << std::endl;

	//Descriptive header:
	ooutput << std::setw(timestepWidth) << "Timestep"
		<< std::setw(timeSecWidth) << "time"
		<< std::setw(fireIntensityWidth) << "Intensity" << std::endl;

	//Units header:
	output << std::setw(timestepWidth) << "Step"
		<< std::setw(timeSecWidth) << "Seconds"
		<< std::setw(fireIntensityWidth) << "kW/m^2" << std::endl;

	//Values:
	for (int i = 0; i < timestep.size(); i++)
	{
		output << std::setw(thickWidth) << timestep[i]//Integer
			<< std::setw(thickWidth) << timeSec[i]//Should be integer.
			<< std::setw(thickWidth) << fireIntensity[i] << std::endl;//Last field don't control the length?
			//<< std::setw(thickWidth) << std::fixed << std::setprecision(2) << fireIntensity[i] << std::endl;
	}

	return output;
}

/** Print the fire history data to an output stream as a set of delimited data rows for each
 * timestep, suitable for data ingestion.
 *
 * @param[in] output The output stream to print to.
 * @param[in] delim The delimiter character.  Defaults to the tab character.
 *
 * @returns The ostream so it can be concatenated to.
 */
std::ostream& BurnupHistory::PrintDelimited(std::ostream& output, const char delim) const
{
	//Print the header:
	output << delim << "Timestep" << delim << "TimeSec" << delim << "FireIntensity" << std::endl;

	//Print the value for each timestep in rows:
	for (int i = 0; i < timestep.size(); i++)
	{
		//The intensity field could be rounded but leave it for accuracy:
		output << delim << timestep[i] << delim << timeSec[i] << delim << fireIntensity[i] << std::endl;
	}
}

//External functions:-------------------------------------------------------------------------------

/** Store the state of a Burnup simulation at the current timestep to a BurnupHistory object for
 * later use.  Sequential calls to this routine will produce a full history of the simulated fire.
 * This is provided as an alternative to saving the history to a file with SaveStateToFile() that
 * make the history programmatically available.
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

/* Overloaded stream print operator for BurnupHistory.
 *
 */
std::ostream& operator<<(std::ostream& output, const BurnupHistory& history)
{
	history.Print(output);
	return output;
}

