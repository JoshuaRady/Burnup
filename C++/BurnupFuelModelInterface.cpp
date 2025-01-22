/***************************************************************************************************
BurnupFuelModelInterface.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 1/17/2025
Reference: Proj. 11 Exp. 22

	This provides an alternate interface to the Burnup wildfire fuel consumption model that uses
fire behavior fuel models, as implemented in the Fireweed Wildfire Code Library.

Licence?????
***************************************************************************************************/

#include <vector>

#include "BurnupFuelModelInterface.h"
#include "BurnupCore.h"
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
 * @param[in] ak			Area influence factor (ak / K_a parameter).
 *              			We modify the original behavior such that a negative value indicates
 *              			that the value of ak / K_a should be calculated according to Albini &
 *              			Reinhardt 1997.
 * @param[in] r0			Minimum value of mixing parameter.
 * @param[in] dr			Max - min value of mixing parameter.
 * @param[in] dT			Time step for integration of burning rates (s). [dT in original Burnup]
 * @param[in] nTimeSteps	Number of time steps to run.
 *
 * @returns A BurnupSim object holding the resulting output from the simulation (and maybe some of the inputs?).
 */
BurnupSim SimulateFM(FuelModel fuelModel,
                     const double duffLoading,// = 0,
                     const double duffMoisture,// = 0,
                     const double tempAirC,
                     const double U,
                     const double fireIntensity,
                     const double t_r,
                     const double ak, const double r0, const double dr, double dT,
                     const int nTimeSteps)
                     //const bool outputHistory = false)///Add?
                     //const bool debug = false);//Add?
{
	const double CtoK = 273;//Burnup's value for 0 C in K.
	
	BurnupSim simData;//Container for simulation data.
	simData.w_o_ij_Initial = fuelModel.w_o_ij;//Store initial fuel loadings.
	int numFuelTypes = fuelModel.numClasses;

	//Check that units of the fuel model are correct:
	if (fuelModel.units != Metric)
	{
		fuelModel.ConvertUnits(Metric);
	}
	
	//Fuel models do not provide names for fuel types so we make some:  Add leading 0s to improve sorting?????
	std::vector<std::string> fuelNames(numFuelTypes);
	
	for (int i = 0; i < numFuelTypes; i++)
	{
		fuelNames[i] = "Fuel " + std::to_string(i);
	}
	
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

	//Other variables that are input only that will be modified...
	double fi = fireIntensity;//Need a copy that can be modified.
	double dtInOut = dT;

	//Set the fuel property inputs missing in the fuel model to default values:
	//In the future these parameters may be added to the FuelModel class or a child thereof.
	//These are default(ish?) values from FOFEM.
	std::vector <double> cheat_ij(numFuelTypes, 2750);//Heat capacity: J/kg K for all fuel types.
	std::vector <double> condry_ij(numFuelTypes, 0.133);//W/m K for all fuel types.
	std::vector <double> tpig_ij(numFuelTypes, 327 + CtoK);//C -> K for intact fuels.
	//FOFEM uses 302 C for punky fuels.
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
	         tpig_ij,//tpig
	         tchar_ij,//tchar,
	         //Outputs by interaction pairs:
	         simData.w_o_kl,//wo
	         simData.xmat_kl,//xmat
	         simData.tign_kl,//tign
	         simData.tout_kl,//tout
	         simData.diam_kl);//diam
	         //outputHistory = false);
	
	//Copy remaining variables to the output object:
	simData.burnoutTime = dtInOut;
	//We could convert negative values to flags.
	
	simData.finalFireIntensity = fireIntensity;//Return the final fire intensity.
	
	//The output by interaction pairs contains useful information but the values by fuel type are
	//most likely to be of primary interest.  Summarize variables by fuel type:
	
	simData.w_o_ij_Final.assign(numFuelTypes, 0);
	
	//Ignition time:
	//Burnup SUMMARY() uses the time the first element in the class ignites.
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

	//Add!!!!!
	//Resort the fuels to their original fuel model order:
	//Since the main outputs are in kl space this may be a bit tricky.

	//Calculate the amount combusted:
	for (int i = 0; i < numFuelTypes; i++)
	{
		simData.combustion_ij[i] = wdry[i] - simData.w_o_ij_Initial[i];
		//simData.combustion_ij[i] = w_o_ij_Final[i] - w_o_ij_Initial[i];//Make sure vectors are in the same order to use this!
	}

	return simData;
}
