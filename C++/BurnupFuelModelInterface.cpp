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
#include "FireweedUnits.h"

/** Perform a fuel consumption simulation using a fuel model and prescribed inputs and return fuel
 * consumption properties.
 *
 * This alternate interface simplifies the inputs and output of Burnup.  ...
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
 *                          (aka R sub M). [dfm in original Burnup]
 *
 * Environmental conditions:
 * @param[in] tempAirC	Ambient air temperature (C).
 * @param[in] U		Mean horizontal windspeed at top of fuelbed [~ at midflame height] (m/s).
 *
 * Igniting fire conditions:
 * @param[in] fireIntensity	...
 * @param[in] t_r			Igniting fire residence time (s). [ti in original Burnup]
 *
 * Simulation conditions and settings:
 * @param[in] ak			Area influence factor (ak / K_a parameter).
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
	const double CtoK = 273;
	
	BurnupSim simData;
	//sim.w_o_ij_Preburn = fm.w_o_ij;//Store initial loading?
	
	//Check that units of the fuel model are correct.  To change it it can be const!
	if (fm.Units != Metric)
	{
		fm.ConvertUnits(Metric);
	}
	
	//Fuel models do not provide names for fuel types so we make some:  Add leading 0s to improve sorting?????
	std::vector<std::string> fuelNames(fuelModel.numClasses);
	
	for (int i = 0; i < fuelModel.numClasses; i++)
	{
		fuelNames[i] = "Fuel " + std::to_string(i);
	}
	
	//Since Burnup may sort the order of fuels we need have to pass in modifiable copies of fuel the
	//properties stored in the fuel model, which is annoying at the least...

	//Fuel model physical properties translated to Burnup terms and units:
	std::vector<double> wdry = fm.w_o_ij;
	std::vector<double> ash = fm.S_T_ij;
	//std::vector<double> htval = fm.h_ij * 1000;//kJ/kg -> J/kg			Need vector math!!!!!
	std::vector<double> htval = fm.h_ij;
	for (double& element : htval)
	{
		element *= 1000;//kJ/kg -> J/kg
	}
	
	std::vector<double> fmois = fm.GetM_f_ij();
	std::vector<double> dendry = fm.rho_p_ij;
	//std::vector<double> sigma = fm.SAV_ij * 100;//cm^2/cm^3 / 1/cm -> m^2/m^3 / 1/m  (1/(m/cm) = cm/m)
	std::vector<double> sigma = fm.SAV_ij;
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
	std::vector <double> cheat_ij(fuelModel.numClasses, 2750);//Heat capacity: J/kg K for all fuel types.
	std::vector <double> condry_ij(fuelModel.numClasses, 0.133);//W/m K for all fuel types.
	std::vector <double> tpig_ij(fuelModel.numClasses, 327 + CtoK);//C -> K for intact fuels.
	std::vector <double> tchar_ij(fuelModel.numClasses, 377 + CtoK);//C -> K for all fuel types.
	
	//Outputs...
	//These don't need to be initialized but currently need to be sized.
// 	BUSim.wo = fm.XXXXX;
// 	BUSim.xmat = fm.XXXXX;
// 	BUSim.tign = fm.XXXXX;
// 	BUSim.tout = fm.XXXXX;
// 	BUSim.diam = fm.XXXXX;
	
	//Perform the simulation:
	Simulate(fireIntensity,//fi,
	         t_r,//ti,
	         U / 60 ,//u: m/s -> m/min
	         fm.delta,//d
	         tempAirC + CtoK,//tpamb C -> K
	         ak, r0, dr, dtInOut,//dt,
	         duffLoading,//df,
	         duffMoisture,//dfm,
	         nTimeSteps,//ntimes,//Move up next to dT?????
	         fm.numClasses;//number
	         
	         fuelNames.//parts
	         wdry, ash, htval, fmois, dendry, sigma,
	         
	         cheat_ij,//cheat
	         condry_ij,//condry,
	         tpig_ij,//tpig
	         tchar_ij,//tchar,
	         //Outputs by interaction pairs:
	         simData.xmat,//xmat
	         simData.tign,//tign
	         simData.tout,//tout
	         simData.w_o_kl,//w_o_out,//wo Note: wo should be moved to the top of the outputs across all interfaces?
	         simData.diam)//diam
	         //outputHistory = false);
	
	//Copy remaining variables to the output object:
	sim.burnoutTime = dtInOut;
	//We could convert negative values to flags.
	
	//sim.finalFireIntensity = fireIntensity;//Return the final fire intensity?
	//sim.w_o_ij_Postburn = //final loadings by fuel type.
	
	
	//Resort the fuels to their original fuel model order:
	//Since the main outputs are in kl space this may be a bit tricky.
	
	
	return sim;
}

