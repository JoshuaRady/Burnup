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


/** Perform a fuel consumption simulation using a fuel model and prescribed inputs and return fuel
 * consumption properties.
 *
 * This alternate interface simplifies the inputs and output of Burnup.  ...
 *
 * @param FuelModel A fuel model ...
 *
 
 * @param[in] U Wind speed [at midflame height] (m/min).
 * @param[in] slopeSteepness				Make sure the units are right!!!!!
 
 * Fuel properties:
 
 *
 * @returns A BUSim object holding the resulting output from the simulation (and maybe some of the inputs?).
 */
BUSim SimulateFM(const FuelModel fuelModel,
                 

                 const double U,// const double slopeSteepness
                 
                 
                 
                 const double fireIntensty,
                 
                 
                 
                 
                 
                 
                 const double ti,
                 //const double u, const double d,
                 const double tpamb,
              const double ak, const double r0, const double dr, double& dt, const double wdf,
              const double dfm,c
              const int ntimes,
              //const int number,
              //std::vector<std::string>& parts,
              //std::vector<double>& wdry, std::vector<double>& ash,
              std::vector<double>& htval, std::vector<double>& fmois, std::vector<double>& dendry,
              std::vector<double>& sigma,
              //std::vector<double>& cheat, std::vector<double>& condry,
              //std::vector<double>& tpig, std::vector<double>& tchar,
              std::vector<double>& xmat,
              std::vector<double>& tign, std::vector<double>& tout, std::vector<double>& wo,
              std::vector<double>& diam,
              const bool outputHistory = false)
                 //const std::vector <double> M_f_ij = {}, const bool debug = false);//Include these?
{
	BUSim output;
	
	
	double fi = fireIntensty;//Need a copy that can be modified.
	
	//Set the fuel property inputs missing in the fuel model to default values:
	std::vector <double> cheat_ij(fuelModel.numClasses, XXXXX);
	std::vector <double> condry_ij(fuelModel.numClasses, XXXXX);
	std::vector <double> tpig_ij(fuelModel.numClasses, XXXXX);
	std::vector <double> tchar_ij(fuelModel.numClasses, XXXXX);
	
	//In the future these parameters may be added to the FuelModel class or a child thereof.
	
	//Since Burnup may sort the order of fuels we need have to pass in modifiable copies of fuel the
	//properties stored in the fuel model, which is annoying at the least...
	std::vector<std::string> parts = XXXXX;
	std::vector<double> wdry =
	std::vector<double> ash = 
	std::vector<double> htval = 
	std::vector<double> fmois = 
	std::vector<double> dendry = 
	std::vector<double> sigma = 
	
	//Other variables that are input only will be modified...
	
	
	BUSim.wo = fm.XXXXX;
	BUSim.xmat = fm.XXXXX;
	BUSim.tign = fm.XXXXX;
	BUSim.tout = fm.XXXXX;
	BUSim.diam = fm.XXXXX;
	
	//Perform the simulation:
	Simulate(fi, const double ti, U,
	         const double d, const double tpamb,
              const double ak, const double r0, const double dr, double& dt, const double wdf,
              const double dfm, const int ntimes, const int number,
              parts,
              wdry, ash,
              htval, fmois, dendry,
              sigma, cheat, condry,
              tpig, tchar, xmat,
              tign, tout, wo,
              diam)//outputHistory = false);
	
	//Copy output to the BUSim object...
	
	
	//Resort the fuels to their fuel model order:
	
	
	//return 1;//Replace with structured output!!!!!
	return output;
}

