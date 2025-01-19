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
 * @param FuelModel A fuel model ...
 *
 
 * @param[in] U Wind speed [at midflame height] (m/min).
 * 				@param[in] slopeSteepness				Make sure the units are right!!!!!
 
 * Fuel properties:
 
 *
 * @returns A BUSim object holding the resulting output from the simulation (and maybe some of the inputs?).
 */
BUSim SimulateFM(//const FuelModel& fuelModel,
                 FuelModel fuelModel,//Pass in by value to the units can be converted if necessary.
                 const double tpamb,
                 const double U,// const double slopeSteepness
                 
                 const double fireIntensty,
                 const double t_r,//ti,
                 
                 //const double u, const double d,
                 //const double tpamb,
                 const double ak, const double r0, const double dr, double dT,
                 const double duffLoading = 0,//wdf,
                 const double duffMoisture = 0,//dfm,
                 const int nTimeSteps,//ntimes,
                 //const int number,
                 //std::vector<std::string>& parts,
                 //std::vector<double>& wdry, std::vector<double>& ash,
                 //std::vector<double>& htval, std::vector<double>& fmois, std::vector<double>& dendry,
                 //std::vector<double>& sigma,
                 //std::vector<double>& cheat, std::vector<double>& condry,
                 //std::vector<double>& tpig, std::vector<double>& tchar,
                 //std::vector<double>& xmat,
                 //std::vector<double>& tign, std::vector<double>& tout, std::vector<double>& wo,
                 //std::vector<double>& diam,
                 const bool outputHistory = false)
                 //const std::vector <double> M_f_ij = {}, const bool debug = false);//Include these?
{
	BUSim sim;//output;
	
	//Check that units of the fuel model are correct.  To change it it can be const!
	if (fm.Units != Metric)
	{
		fm.ConvertUnits(Metric);
	}
	
	//Fuel models do not provide names for fuel types so we make some:  Add leading 0s to improve sorting?????
	std::vector<std::string> parts(fuelModel.numClasses);//fuelNames?????
	
	for (int i = 0; i < *number; i++)
	{
		parts[i] = "Fuel " + std::to_string(i);
	}
	
	//Since Burnup may sort the order of fuels we need have to pass in modifiable copies of fuel the
	//properties stored in the fuel model, which is annoying at the least...

	//Fuel model physical properties translated to Burnup terms and units:
	std::vector<double> wdry = fm.w_o_ij;
	std::vector<double> ash = fm.S_T_ij;
	std::vector<double> htval = fm.h_ij * 1000;//kJ/kg -> J/kg			Need vector math!!!!!
	std::vector<double> fmois = fm.GetM_f_ij();
	std::vector<double> dendry = fm.rho_p_ij;
	std::vector<double> sigma = fm.SAV_ij * 100;//cm^2/cm^3 / 1/cm -> m^2/m^3 / 1/m  (1/(m/cm) = cm/m)
	
	//Other variables that are input only that will be modified...
	double fi = fireIntensty;//Need a copy that can be modified.
	double dtInOut = dT;
	
	
	//Set the fuel property inputs missing in the fuel model to default values:
	//In the future these parameters may be added to the FuelModel class or a child thereof.
	//These are default(ish?) values from FOFEM.
	std::vector <double> cheat_ij(fuelModel.numClasses, 2750);//Heat capacity: J/kg K for all fuel types.
	std::vector <double> condry_ij(fuelModel.numClasses, 0.133);//W/m K for all fuel types.
	std::vector <double> tpig_ij(fuelModel.numClasses, 327 + 273);//C -> K for intact fuels.
	std::vector <double> tchar_ij(fuelModel.numClasses, 377 + 273);//C -> K for all fuel types.
	
	//Outputs...
	//These don't need to be initialized but currently need to be sized.
// 	BUSim.wo = fm.XXXXX;
// 	BUSim.xmat = fm.XXXXX;
// 	BUSim.tign = fm.XXXXX;
// 	BUSim.tout = fm.XXXXX;
// 	BUSim.diam = fm.XXXXX;
	
	//Perform the simulation:
	Simulate(fireIntensty,//fi,
	         t_r,//ti,
	         U,//u
	         fm.delta,//d
	         const double tpamb,
	         ak, r0, dr, dtInOut,//dt,
	         duffLoading,//df,
	         duffMoisture,//dfm,
	         nTimeSteps,//ntimes,//Move up next to dT?????
	         fm.numClasses;//number
	         
	         parts,
	         wdry, ash, htval, fmois, dendry, sigma,
	         
	         cheat_ij,//cheat
	         condry_ij,//condry,
	         tpig_ij,//tpig
	         tchar_ij,//tchar,
	         //Outputs:
	         sim.xmat,//xmat
	         sim.tign,//tign
	         sim.tout,//tout
	         sim.w_o_kl,//w_o_out,//wo Note: wo should be moved to the top of the outputs across all interfaces?
	         sim.diam)//diam
	         //outputHistory = false);
	
	//Copy output to the BUSim object...
	
	
	//Resort the fuels to their fuel model order:
	
	
	//return 1;//Replace with structured output!!!!!
	return output;
}

