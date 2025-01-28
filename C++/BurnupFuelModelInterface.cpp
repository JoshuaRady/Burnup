/***************************************************************************************************
BurnupFuelModelInterface.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 1/17/2025
Reference: Proj. 11 Exp. 22

	This provides an alternate interface to the Burnup wildfire fuel consumption model that uses
fire behavior fuel models, as implemented in the Fireweed Wildfire Code Library.

References:
Albini, F.A., Brown, J.K., Reinhardt, E.D. and Ottmar, R.D.
Calibration of a large fuel burnout model.
International Journal of Wildland Fire 5(3): 173-192, 1995.

	This paper provides insight into setting parameters for Burnup.

Licence?????
***************************************************************************************************/

#include <vector>
#include <iostream>
#include <iomanip>

#include "BurnupFuelModelInterface.h"
#include "BurnupCore.h"
#include "FireweedMessaging.h"
#include "FireweedStringUtils.h"//For PrintVector().
#include "FireweedUnits.h"

/** Perform a fuel consumption simulation using a fuel model and prescribed inputs and return fuel
 * consumption properties.
 *
 * This alternate interface allows a fuel model to be used as input to Burnup.  This reduces the
 * number f inputs and simplifies outputs by returning all of them in a single object rather than
 * using multiple return arguments which need to be initialized by the calling code.  The parameter
 * names have been made more descriptive than the standard interface.  The code handles conversion
 * of variables that differ in units between the standard fuel models and Burnup.
 *
 * Fuel models do not contain all the fuel properties that Burnup needs.  For now the heat capacity
 * (cheat), thermal conductivity (condry), ignition temperature (tpig), and char temperature (tchar)
 * for all fuels are set to default values.  In the future these parameters may be added to the
 * FuelModel class or a child thereof.
 *
 * Fuel properties:
 * @param FuelModel A fuel model to calculate fuel consumption for.
 *                  The fuel moisture, M_f_ij, must be included in the FuelModel object.
 *                  FuelModel is not const so we can the units can be converted if necessary.
 *                  The post-fire loading could be returned but updating w_o_ij but given the many
 *                  other outputs it doesn't seem useful.
 * 
 * Duff conditions:
 * @param[in] duffLoading	Duff loading (kg/m^2, aka W sub d). [wdf in original Burnup]
 * @param[in] duffMoisture	Ratio of moisture mass to dry organic mass / duff fractional moisture
 *                        	aka R sub M). [dfm in original Burnup]
 *
 * Environmental conditions:
 * @param[in] tempAirC	Ambient air temperature (C).
 * @param[in] U			Mean horizontal windspeed at top of fuelbed [~ at midflame height] (m/s).
 *
 * Igniting fire conditions:
 * @param[in] fireIntensity	Igniting fire intensity (site avg) (kW / m^2).
 * @param[in] t_r			Igniting fire residence time (s). [ti in original Burnup]
 *
 * Simulation conditions and settings:
 * @param[in] dT			Time step for integration of burning rates (s). [dT in original Burnup]
 * @param[in] nTimeSteps	Number of time steps to run.
 * @param[in] burnupFormat	How should the output be returned?  The default (false) is to return the
 *                        	fuel level output in the same units and fuel order as the fuel model
 *                        	input.  If true the output will use the original Burnup units and the
 *                        	fuel level output may be reordered.  In either case the outputs at the
 *                        	fuel pair level will be raw Burnup output (at least for now).
 *
 * @note The order for the optional parameters is tricky.  We try to put them in the order with the
 * most likely to be supplied first.
 *
 * Optional fuel properties:
 * @param[in] tpig_ij		Ignition temperature, C.  Defaults to 327 C for all fuels.
 *
 * Optional conditions and settings:
 * @param[in] ak			Area influence factor (ak / K_a parameter).
 *              			We modify the original behavior such that a negative value indicates
 *              			that the value of ak / K_a should be calculated according to Albini &
 *              			Reinhardt 1997.  By default calculate the value.
 * @param[in] r0			Minimum value of mixing parameter.
 *              			Default value of 1.83 from Albini, Brown, Reinhardt, and Ottmar 1995.
 * @param[in] dr			Max - min value of mixing parameter.
 *              			Default value of 0.40 from Albini, Brown, Reinhardt, and Ottmar 1995.
 *
 * @returns A BurnupSim object holding the resulting output from the simulation (and maybe some of the inputs?).
 *
 * @note The setting values can be hard to determine.  We provide default values for these based on
 * on some papers but more could be done to inform value selection.
 * @note Setting duffLoading = 0 should eliminate the effect of duff on the simulation.  Duff
 * moisture should not matter with no loading.  1.0 seems like a good placeholder value in this
 * case.
 */
BurnupSim SimulateFM(FuelModel fuelModel,
                     const double duffLoading,//Default to 0?
                     const double duffMoisture,//Default to 1.0?
                     const double tempAirC, const double U,
                     const double fireIntensity, const double t_r,
                     const double dT, const int nTimeSteps,
                     const bool burnupFormat,
                     const std::vector <double> tpig_ij,
                     const double ak, const double r0, const double dr)
                     //const bool outputHistory = false)///Add?
                     //const bool debug = false);//Add?
{
	const double CtoK = 273;//Burnup's value for 0 C in K.
	
	BurnupSim simData;//Container for simulation data.
	
	int numFuelTypes = fuelModel.numClasses;

	//Check that units of the fuel model are metric as needed:
	if (fuelModel.units != Metric)
	{
		fuelModel.ConvertUnits(Metric);
	}
	
	//Fuel models do not provide names for fuel types so we make some:
	//Add leading 0s to improve sorting?????
	std::vector<std::string> fuelNames(numFuelTypes);
	
	for (int i = 0; i < numFuelTypes; i++)
	{
		fuelNames[i] = "Fuel " + std::to_string(i + 1);
	}
	std::vector<std::string> fuelNamesInitial = fuelNames;//fuelNamesFM??????
	
	//Since Burnup may sort the order of fuels we need to pass in modifiable copies of fuel
	//properties stored in the fuel model:

	//Fuel model physical properties translated to Burnup terms and units:
	std::vector<double> wdry = fuelModel.w_o_ij;
	std::vector<double> ash = fuelModel.S_T_ij;
	std::vector<double> htval = fuelModel.h_ij;
	for (double& element : htval)
	{
		element *= 1000;//kJ/kg -> J/kg
	}

	std::vector<double> fmois = fuelModel.GetM_f_ij();
	std::vector<double> dendry = fuelModel.rho_p_ij;
	std::vector<double> sigma = fuelModel.SAV_ij;
	for (double& theSAV : sigma)
	{
		theSAV *= 100;//cm^2/cm^3 = 1/cm -> m^2/m^3 = 1/m  (1/(m/cm) = cm/m = 100)
	}

	//Make copies of other variables that are input only that will be modified by Burnup:
	double fi = fireIntensity;
	double dtInOut = dT;

	//Set the fuel property inputs missing in the fuel model to default values:
	//In the future these parameters may be added to the FuelModel class or a child thereof.
	//These are default(ish?) values from FOFEM.
	std::vector <double> cheat_ij(numFuelTypes, 2750);//Heat capacity: J/kg K for all fuel types.
	std::vector <double> condry_ij(numFuelTypes, 0.133);//Conductivity: W/m K for all fuel types.
	
	//std::vector <double> tpig_ij(numFuelTypes, 327 + CtoK);//Ignition Temp: C -> K for intact fuels.
	std::vector <double> tpig_ij_K(numFuelTypes);//Modifiable copy in Kelvins.
	if (tpig_ij.empty())
	{
		tpig_ij_K.assign(numFuelTypes, 327 + CtoK);//Ignition Temp: C -> K for intact fuels.
		//FOFEM uses 302 C for 'punky' logs.
	}
	else
	{
		//We could allow a single value be passed in as a value for all members.
		if (tpig_ij.size() != numFuelTypes)
		{
			Stop("tpig_ij does not match the number of fuel types.");
		}

		//for (int i = 0; i < numFuelTypes; i++)
		//for (int k = 0; k < numFuelTypes; k++)//k used here in the fuel model sense.
		for (int kFM = 0; kFM < numFuelTypes; kFM++)//k used here in the fuel model sense.
		{
			tpig_ij_K[kFM] = tpig_ij[kFM] + CtoK;
		}
	}
	
	std::vector <double> tchar_ij(numFuelTypes, 377 + CtoK);//C -> K for all fuel types.

	//Perform the simulation:
	Simulate(fi,//fi,
	         t_r,//ti,
	         U / 60 ,//u: m/min -> m/s
	         fuelModel.delta,//d
	         tempAirC + CtoK,//tpamb C -> K
	         ak, r0, dr, dtInOut,//dt,
	         duffLoading,//df,
	         duffMoisture,//dfm,
	         nTimeSteps,//ntimes,//Move up next to dT?????
	         numFuelTypes,//number
	         //Fuel properties:
	         fuelNames,//parts
	         wdry, ash, htval, fmois, dendry, sigma,
	         cheat_ij,//cheat
	         condry_ij,//condry,
	         tpig_ij_K,//tpig
	         tchar_ij,//tchar,
	         //Outputs by interaction pairs:
	         simData.w_o_kl,//wo
	         simData.xmat_kl,//xmat
	         simData.tign_kl,//tign
	         simData.tout_kl,//tout
	         simData.diam_kl);//diam
	         //outputHistory = false);

	//Copy output variables [not directly modified by Burnup] to the output object:
	simData.burnoutTime = dtInOut;
	//We could convert negative values to flags.

	simData.finalFireIntensity = fireIntensity;//Return the final fire intensity.

	simData.klFuelNames = fuelNames;//Alway store the potentially reordered names.


	//The outputs by interaction pairs contain useful information but the values by fuel type are
	//most likely to be of primary interest.  Summarize variables by fuel type:
	simData.w_o_ij_Final.assign(numFuelTypes, 0);
	
	//Ignition times:
	//Burnup SUMMARY() uses the time the first element in the class ignites so we do too.
	simData.tign_ij.assign(numFuelTypes, 0);
	
	//Burnout time:
	//For some pairs the burnout time may not have been computed and the value will be rindef.  In
	//the original code SUMMARY() just takes the max including these values.  That doesen't make a
	//lot of sense.  We calculate the minimum and maximum for each fuel type:
	simData.tout_ij_Min.assign(numFuelTypes, 0);
	simData.tout_ij_Max.assign(numFuelTypes, 0);
	
	for (int k = 1; k <= numFuelTypes; k++)
	{
		int k0 = k - 1;
		
		for (int l = 0; l <= k; l++)//l in kl space, 0 based
		{
			int kl = Loc(k, l);
			simData.w_o_ij_Final[k0] += simData.w_o_kl[kl];
			
			if (l == 0)
			{
				simData.tign_ij[k0] = simData.tign_kl[kl];
				simData.tout_ij_Min[k0] = simData.tout_kl[kl];
			}
			else
			{
				simData.tign_ij[k0] = std::min(simData.tign_ij[k0], simData.tign_kl[kl]);
				//This should ignore rindef:
				simData.tout_ij_Min[k0] = std::min(simData.tout_ij_Min[k0], simData.tout_kl[kl]);
			}
			
			simData.tout_ij_Max[k0] = std::max(simData.tout_ij_Max[k0], simData.tout_kl[kl]);
		}

		//Make sure no burnout time exceeds that of the actual fire:
		//Estimates or rindef values may lead to this.
		simData.tout_ij_Min[k0] = std::min(simData.tout_ij_Min[k0], simData.burnoutTime);
		simData.tout_ij_Max[k0] = std::min(simData.tout_ij_Max[k0], simData.burnoutTime);
	}

	//This needs to be tested for cases where some fuels don't ignite!!!!

	if (burnupFormat)
	{
		simData.fuelModelFormat = false;

		//Store the fuel level inputs in their reordered (and re-united) forms:
		simData.SAV_ij = sigma;//SAVs.
		simData.w_o_ij_Initial = wdry;//Initial fuel loadings.
		simData.M_f_ij = fmois;//Fuel moisture/

		simData.fuelNames = fuelNames;
	}
	else
	{
		simData.fuelModelFormat = true;

		//Store the fuel level inputs in their original fuel model order and metric units:
		simData.SAV_ij = fuelModel.SAV_ij;//SAVs
		simData.w_o_ij_Initial = fuelModel.w_o_ij;//Initial fuel loadings.
		simData.M_f_ij = fuelModel.GetM_f_ij();//Fuel moisture.

		simData.fuelNames = fuelNamesInitial;

		std::vector<int> fuelOrder(numFuelTypes);

		if (fuelNamesInitial == fuelNames)
		{
			//fuelOrder = ?????
			simData.klFuelsReordered = false;
			//There is no need tp resort.
		}
		else
		{
			simData.klFuelsReordered = true;

			//Determine how the fuels were reordered:
			//Burnup stores the sort order internally as key[].  We could pass that out to eliminate the
			//need for this code.
			for (int m = 0; m < numFuelTypes; m++)
			{
				for (int n = 0; n < numFuelTypes; n++)
				{
					if (fuelNamesInitial[m] == fuelNames[n])
					{
						//This gives the the postions each original index was move to:
						//fuelOrder[m] = n;
						//What we want to is were to move indexes to get the original order back:
						fuelOrder[n] = m;
						break;
					}
				}
			}

			//Reorder the fuel level outputs to the original fuel model order:
			Reorder(simData.w_o_ij_Final, fuelOrder);
			Reorder(simData.tign_ij, fuelOrder);
			Reorder(simData.tout_ij_Min, fuelOrder);
			Reorder(simData.tout_ij_Max, fuelOrder);

			//I don't think any output units need to be converted?????

			//The outputs by pairs are not easily reordered so we leave them for now.
		}
	}

	//Calculate the amount combusted at the end so we don't have to reorder it:
	simData.combustion_ij.assign(numFuelTypes, 0);
	for (int i = 0; i < numFuelTypes; i++)
	{
		simData.combustion_ij[i] = simData.w_o_ij_Initial[i] - simData.w_o_ij_Final[i];

		//Add checking for values that are below 0?
	}

	//Need to add handling for history data!!!!!

	return simData;
}

/** Reorder a vector.
 * 
 * @param vec	The vector to reorder.
 * @param order	The order (0 based) that the current indexes should be moved to.
 *
 * @returns Nothing.  We could return the vector rather than modifying it in place.
 */
void Reorder(std::vector<double>& vec, const std::vector<int> order)
{
	//Check the order is valid:
	if (vec.size() != order.size())
	{
		Stop("Reorder(): The order is not the same length as the vector.");
	}
	//Check order values are valid...

	//std::vector<double> temp(vec.size())
	std::vector<double> copy = vec;

	for (int i = 0; i < vec.size(); i++)
	{
		vec[order[i]] = copy[i];
	}
}

//BurnupSim Functions:------------------------------------------------------------------------------

/** Print the Burnup simulation data to an output stream.
 *
 * @param output The output stream to print to.
 *
 * @returns The ostream so it can be conatinated to.
 *
 * @note This is a draft and only some of the top level data members are currently printed.
 */
std::ostream& BurnupSim::Print(std::ostream& output) const
{
	output << "Burnup simulation:" << std::endl;

	if (fuelModelFormat)
	{
		output << "Data is in fuel moddel format." << std::endl;
	}
	else
	{
		output << "Data is in Burnup format." << std::endl;
	}

	if (burnoutTime == -1.0)
	{
		output << "Igniting fire cannot ignite fuel." << std::endl;
	}
	else if (burnoutTime == -2.0)
	{
		output << "Igniting fire cannot dry fuel." << std::endl;
	}
	else
	{
		output << "The fire burnt out after " << burnoutTime << " seconds." << std::endl;
		output << "The final fire intensity was " << finalFireIntensity << " (kW / m^2)." << std::endl;

		//Print one by one:
// 		output << "Fuel Names: ";
// 		//PrintVector(output, fuelNames);
// 		for (int i = 0; i < fuelNames.size() - 1; i++)
// 		{
// 			output << fuelNames[i] << ", ";
// 		}
// 		output << fuelNames[fuelNames.size() - 1] << std::endl;
// 
// 		output << "SAV_ij: ";
// 		PrintVector(output, SAV_ij);
// 
// 		output << "M_f_ij: ";
// 		PrintVector(output, M_f_ij);
// 
// 		output << "w_o_ij_Initial: ";
// 		PrintVector(output, w_o_ij_Initial);
// 
// 		output << "w_o_ij_Final: ";
// 		PrintVector(output, w_o_ij_Final);
// 
// 		output << "combustion_ij: ";
// 		//output << "Fuel combusted: ";
// 		PrintVector(output, combustion_ij);
// 
// 		output << "tign_ij: ";
// 		PrintVector(output, tign_ij);
// 
// 		output << "tout_ij_Min: ";
// 		PrintVector(output, tout_ij_Min);
// 
// 		output << "tout_ij_Max: ";
// 		PrintVector(output, tout_ij_Max);

		//Print the the fuel level outputs in a table for easy interpretation:

		//Widths similar to standard Burnup output table formating elsewhere:
		//const int nameWidth = 7;
		//const int w_oIWidth = 11;
		//const int w_oFWidth = 12;
		//const int tignWidth = 16;
		//const int toutMinWidth = 11;
		//const int toutMaxWidth = 11;
		//const int m_fWidth = 10;
		//const int savWidth = 8;

		const int nameWidth = 7;
		const int w_oIWidth = 15;
		const int w_oFWidth = 13;
		const int tignWidth = 16;
		const int toutMinWidth = 12;
		const int toutMaxWidth = 12;
		const int m_fWidth = 10;
		int savWidth;
		std::string savUnits;
		
		if (fuelModelFormat)
		{
			savWidth = 8;//9;
			savUnits = "cm^2/cm^3";
		}
		else
		{
			savWidth = 8;
			savUnits = "m^2/m^3";
		}

		//Member name header:
		std::cout << std::setw(nameWidth) << "Name"
			<< std::setw(w_oIWidth) << "w_o_ij_Initial"
			<< std::setw(w_oFWidth) << "w_o_ij_Final"
			<< std::setw(tignWidth) << "tign_ij"
			<< std::setw(toutMinWidth) << "tout_ij_Min"
			<< std::setw(toutMaxWidth) << "tout_ij_Max"
			<< std::setw(m_fWidth) << "M_f_ij"
			<< std::setw(savWidth) << "SAV_ij" << std::endl;

		//Units header:
		std::cout << std::setw(nameWidth) << " "
			<< std::setw(w_oIWidth) << "kg/m^2"
			<< std::setw(w_oFWidth) << "kg/m^2"
			<< std::setw(tignWidth) << "seconds"
			<< std::setw(toutMinWidth) << "seconds"
			<< std::setw(toutMaxWidth) << "seconds"
			<< std::setw(m_fWidth) << "fraction"
			<< std::setw(savWidth) << savUnits << std::endl;

		//Descriptive header:
		std::cout << std::setw(nameWidth) << " "
			<< std::setw(w_oIWidth) << "Preburn_wo"
			<< std::setw(w_oFWidth) << "Postburn_wo"
			<< std::setw(tignWidth) << "IgnitionTimeMin"
			<< std::setw(toutMinWidth) << "BurnoutMin"
			<< std::setw(toutMaxWidth) << "BurnoutMax"
			<< std::setw(m_fWidth) << "MoistFrac"
			<< std::setw(savWidth) << "SAV" << std::endl;

		//Values:
		std::streamsize thePrecision = std::cout.precision();
		for (int i = 0; i < SAV_ij.size(); i++)
		{
			std::cout << std::setw(nameWidth) << fuelNames[i]
				<< std::setw(w_oIWidth) << std::fixed << std::setprecision(5) << w_o_ij_Initial[i]
				<< std::setw(w_oFWidth) << std::fixed << std::setprecision(5) << w_o_ij_Final[i]
				<< std::setw(tignWidth) << std::fixed << std::setprecision(0) << tign_ij[i]
				<< std::setw(toutMinWidth) << std::fixed << std::setprecision(0) << tout_ij_Min[i]
				<< std::setw(toutMaxWidth) << std::fixed << std::setprecision(0) << tout_ij_Max[i]
				<< std::setw(m_fWidth) << std::fixed << std::setprecision(2) << M_f_ij[i]
				<< std::setw(savWidth) << std::fixed << std::setprecision(2) << SAV_ij[i] << std::endl;
		}
		std::cout.precision(thePrecision);//Restore the previous setting.
	}

	return output;
}

/* Overloaded stream print operator for BurnupSim.
 *
 */
std::ostream& operator<<(std::ostream& output, const BurnupSim& buSim)
{
	buSim.Print(output);
	return output;
}
