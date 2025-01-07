/***************************************************************************************************
Burnup.cpp
Burnup Wildfire Fuel Consumption Model

Original Fortan code by: Frank A. Albini
Edited, modified, and ported to C++ by: Joshua M. Rady
Woodwell Climate Research Center
Started: 12/30/2024
Reference: Proj. 11 Exp. 22

This is an reimplementation of the Burnup wildfire fuel consumption model in C++.

...

***************************************************************************************************/

#include <algorithm>//For max(), min(), fill().
#include <cerrno>
#include <cmath>//<math.h>//For pow(), sqrt(), abs().
#include <cstring>//For strerror().
#include <fstream>
#include <iostream>//Or just <ostream>?
#include <vector>

#include "Burnup.h"
#include "FireweedMessaging.h"

/* Program level dimensional constants:
In the original code these parameters were passed into all routines that need them.  Here we
change to module level scope.  This allows us to reduce the number of arguments to routines.
While currently fixed, these can probably be made dynamic to be set a initialization.*/
//Update!!!!!

/*The maximum number of fuel components or types.  This is used to build fixed size data
structures.  The number of elements may be less than this so some indexes may be empty.
The original program fixes this arbitrarily at 10 fuel components.*/
const int maxno = 12;
//The maximum number of non-zero entries in the triangular matrix of fuel interaction pairs:
//Add one to one dimension for the  'no companion' interaction element.
const int maxkl = maxno * (maxno + 1) / 2 + maxno;
//The maximum dimension of historical sequences (primarily for qdot):
const int mxstep = 20;

//Physical constants:
//In the original code these were defined as variables and shared through subroutine arguments.
const double ch2o = 4186.0;//Specific heat capacity of water, J / kg K
const double tpdry = 353.0;//Temperature (all components) start drying (K)

//Empty / undefined value constant:
const double rindef = 1.0e+30;//In the original code this defined both in START() and STEP().

//Names:
const char noCmpStr[] = "no companion";//The name for no companion pairs.
                                       //In the original code this was declared as a variable 'none'
                                       //in Summary().  'none' is a Fortran keyword so it was renamed.

//Globals:
bool SaveHistory = false;//Should fire history be output to file?


//InteractiveUI()

//Simulate()
/** Perform a simulation with prescribed inputs and return fuel consumption properties.
 * The main fire properties are returned as output parameters.  Optionally additional a detailed
 * fire history can be output to file.
 *
 * Igniting fire and environmental data:
 * @param fi	Current fire intensity (site avg), kW / sq m
 * The value passed in for fi is the fire front intensity.  [Actually this is probably not quite right!!!!!  It is probably not the FRONT intensity...]
 * The variable is later reused and
 * updated by FIRINT().  It is passed on to other routines that use but do not change it.
 * These two uses could be separated.  The value returned is the final intensity, which might
 * be of use.  A history would be more valuable.
 *
 * @param ti		Igniting fire residence time (s).
 * @param u			Mean horizontal windspeed at top of fuelbed (m/s).
 * @param d			Fuelbed depth (m).
 * @param tpamb		Ambient temperature (K).
 * 
 * Internal and control variables:
 * @param ak		Area influence factor (ak / K_a parameter).
 * @param r0		Minimum value of mixing parameter.
 * @param dr		Max - min value of mixing parameter.
 * @param dt		Time step for integration of burning rates (s).
 *          		On completion contains the time the fire went out.
 *          		A value of -1 indicates the fuel did not ignite.  A value
 *          		of -2 indicates the fuel did not complete drying.  In such
 *          		cases most of remaining return variables will not be
 *          		meaningful.
 * 
 * Considering removing these two.  See below:
 * @param wdf		Duff loading (kg/m^2, aka W sub d).
 * @param dfm		Ratio of moisture mass to dry organic mass /
 *           		duff fractional moisture (aka R sub M).
 * @param ntimes	Number of time steps to run.  Move down?
 * @param number	The number of fuel classes.
 *
 * Fuel component property arrays:  The values will not change but they may be reordered.
 * Returning the reordered arrays may be overkill.  The revised order might be sufficient.
 * However, setting these to inout allows the values to be reordered internally by SORTER(),
 * which eliminates the need for parallel local variables.
 * @param parts		Fuel component names / labels. [maxno]
 * @param wdry		Ovendry mass loading, kg/sq m. [maxno]
 * @param ash		Mineral content, fraction dry mass. [maxno]
 * @param htval		Low heat of combustion (AKA heat content), J / kg. [maxno]
 * @param fmois		Moisture fraction of component. [maxno]
 * @param dendry	Ovendry mass density, kg / cu m. [maxno]
 * @param sigma		Surface to volume ratio, 1 / m. [maxno]
 * @param cheat		Specific heat capacity, (J / K) / kg dry mass. [maxno]
 * @param condry	Thermal conductivity, W / m K, ovendry. [maxno]
 * @param tpig		Ignition temperature, K. [maxno]
 * @param tchar		Char temperature, K. [maxno]
 *
 * Calculated outputs:
 * The following are the main variables output by SUMMARY(): [name], fr, ti, to, wd, di
 *
 * Use the first argument with maxkl length to calculate its value.  This can then be used
 * by subsequent variables.  I don't like this approach much but a better alternative has
 * not be determined.
 * wo should be moved to the front in any case because it is the most valuable output.  This
 * would have the advantage of making the size() shorthand shorter.
 * JMR_Note: No longer in argument order!!!!!
	! real*4, intent(out) :: wo(number * (number + 1) / 2 + number)
 * @param wo		Current ovendry loading for the larger of each component pair, kg/sq m. [maxkl]
 * @param xmat		Table of influence fractions between components. [maxkl]
 * @param tign		Ignition time for the larger of each fuel component pair. [maxkl]
 * @param tout		Burnout time of larger component of pairs. [maxkl]
 * @param diam		Current diameter of the larger of each fuel component pair, m. [maxkl]
 *
 * @note Currently it is expected that the calculated output vectors be supplied in the correct size.
 * We could move the sizing into the function allowing empty vectors to be passed in.
 *
 * @returns On return many arguments may be sorted, updated, or returned. [MORE!!!!!]
 *
 * Settings:
 * @param outputHistory	Should fire history be saved?  Defaults to false.
 *

! History: Added as an programatic alternative entry point to the original interactive program.
 */
void Simulate(double& fi, const double ti, const double u, const double d, const double tpamb,
              const double ak, const double r0, const double dr, double& dt, const double wdf,
              const double dfm, const int ntimes, const int number,
              std::vector<std::string>& parts, std::vector<double>& wdry, std::vector<double>& ash,
              std::vector<double>& htval, std::vector<double>& fmois, std::vector<double>& dendry,
              std::vector<double>& sigma, std::vector<double>& cheat, std::vector<double>& condry,
              std::vector<double>& tpig, std::vector<double>& tchar, std::vector<double>& xmat,
              std::vector<double>& tign, std::vector<double>& tout, std::vector<double>& wo,
              std::vector<double>& diam, const bool outputHistory)
{
	//Arrays:
	std::vector<double> alfa(number);		//Dry thermal diffusivity of component, sq m / s. [maxno]
	std::vector<double> flit(number);		//Fraction of each component currently alight. [maxno]
	std::vector<double> fout(number);		//Fraction of each component currently gone out. [maxno]
	std::vector<double> work(number);		//Workspace array. [maxno]
	std::vector<std::vector<double>> elam(number, std::vector<double>(number));//Interaction matrix. [maxno, maxno]
	std::vector<double> alone(number);		//Non-interacting fraction for each fuel class. [maxno]
	std::vector<double> area(number);		//Fraction of site area expected to be covered at
	                                 		//least once by initial planform area of ea size. [maxno]
	std::vector<double> fint(number);		//Corrected local fire intensity for each fuel type. [maxno]
	std::vector<double> tdry(wo.size());	//Time of drying start of the larger of each fuel component pair. [maxkl]
	std::vector<double> wodot(wo.size());	//Dry loading loss rate for larger of pair. [maxkl]
	std::vector<double> ddot(wo.size());	//Diameter reduction rate, larger of pair, m / s. [maxkl]
	std::vector<double> qcum(wo.size());	//Cumulative heat input to larger of pair, J / sq m. [maxkl]
	std::vector<double> tcum(wo.size());	//Cumulative temp integral for qcum (drying). [maxkl]
	std::vector<double> acum(wo.size());	//Heat pulse area for historical rate averaging. [maxkl]
	std::vector<std::vector<double>> qdot(wo.size(), std::vector<double>(mxstep));//History (post ignite) of heat transfer rate
	                                    	//to the larger of each component pair, W / sq m. [maxkl, mxstep]
	std::vector<int> key(number);			//Ordered index list. [maxno]
	std::vector<std::string> list(number);		//Intermediary for reordering parts name array. [maxno]
	                                 		//Probably not needed here.  See notes in ARRAYS().

	//Scalars in order of appearance:
	double dfi; 			//Duff fire intensity (aka I sub d) for DUFBRN().
	double tdf; 			//Burning duration (aka t sub d) for DUFBRN().
	int now;			//Index marking the current time step.
	double tis;				//Current time (ti + number of time steps * dt).
	int ncalls;		//Counter of calls to START().
	double fid;				//Fire intensity due to duff burning.

	//In the original code fmin was a local treated as a constant.  Passing it in might be good:
	const double fimin = 0.1;//Fire intensity (kW / sq m) at which fire goes out.

	//There are a large number of locals in this routine that are not explictly initialized.
	//Most are initialized in ARRAYS() and START().  Testing was done to confrim that explicit
	//initialization was not needed for the remaining locals.

	//Set SaveHistory:
// 	if (present(outputHistory)) then
// 		SaveHistory = outputHistory
// 	else
// 		SaveHistory = .false.
// 	end if
	SaveHistory = outputHistory;

	//Adding size checking on incoming arrays would be good here.

	//Sort the fuel components and calculate the interaction matrix...
	ARRAYS(wdry, ash, dendry, fmois, sigma, htval, cheat, condry, tpig, tchar, diam, key, work, ak,
	       elam, alone, xmat, wo, parts, list, area);

	//Record the state before the start of the simulation.  This need to be done after ARRAYS()
	//because parts, wo, and diam may get reordered.  now and tis are not initialized yet and we
	//set the time explicitly.  
	//SaveStateToFile(0, 0.0, number, parts, wo, diam, fi);
	//The first simulated time point starts at the end of the Igniting fire residence time.  The
	//fire intensity is a constant value for this period.  It would make sense to make another
	//record of fire intensity at the end of the residence time.  However, this would lead to
	//two values at time point 1.  Likewise, this would result in two values for fuel loadings.
	//Both fire intensity and fine fuel loads can drop significantly at the first time step.

	//The original code calls DUFBRN() here.  I'm leaving this here while getting the code
	//running but it would probably better to pass the output of DUFBRN() in instead.
	DUFBRN(wdf, dfm, dfi, tdf);

	//Initialize variables and data structures:
	now = 1;
	tis = ti;
	START(tis, now, wo, alfa, dendry, fmois, cheat, condry, diam, tpig, tchar, xmat, tpamb, fi,
	      flit, fout, tdry, tign, tout, qcum, tcum, acum, qdot, ddot, wodot, work, u, d, r0, dr,
	      ncalls);

	//if (tign(1) .lt. 0.0) then //Fuels failed to ignite.
	if (tign[0] < 0.0)//Fuels failed to ignite.
	{
		//We could use dt = 0 to indicate the condition but could be confused as an actual
		//value.  A negative value is clearly not valid and the value can be used to provide
		//more information.
		dt = tign[0];
		return;
	}

	//If the duff burns longer than the passing fire front then have it's intensity
	//contribute to the post front fire environment, otherwise ignore it:
	if (tis < tdf)
	{
		fid = dfi;
	}
	else
	{
		fid = 0.0;
	}

	//Calculate the initial fire intensity:
	FIRINT(wodot, ash, htval, number, area, fint, fi);

	//Record the state after START() and the first call to FIRINT(): Make optional!!!!!
	//SaveStateToFile(now, tis, number, parts, wo, diam, fi);

	//If the fire intensity is above the extinguishing threshold calculate combustion until
	//the fire goes out or the number of timesteps is reached:
	if (fi > fimin)
	{
		while (now < ntimes)
		{
			STEP(dt, now, wo, alfa, dendry, fmois, cheat, condry, diam, tpig, tchar, xmat, tpamb,
			     fi, flit, fout, tdry, tign, tout, qcum, tcum, acum, qdot, ddot, wodot, work, u, d,
			     r0, dr, ncalls, tis, fint, fid);

			//Update time trackers:
			now += 1;
			tis = tis + dt;
			//JMR_NOTE: It is a bit strange that START() and STEP() both start at the same
			//time step.  Think on this a bit!!!!!

			//Get the duff contribution while it remains burning:
			if (tis < tdf)
			{
				fid = dfi;
			}
			else
			{
				fid = 0.0;
			}

			//Calculate the fire intensity at this time step:
			FIRINT(wodot, ash, htval, number, area, fint, fi);

			//Save the state at each timestep: Make optional!!!!!
			//SaveStateToFile(now, tis, number, parts, wo, diam, fi);

			if (fi <= fimin)
			{
				//exit
				break;
			}
		}
	}//(fi > fimin)

	//Return the time when the fire dropped below fimin as the time the fire "went out".
	//We return the value in dt.  This may be changed in the future.
	dt = tis;
	//We could also return the number of timesteps completed in ntimes but that doesn't add much.
}

//SimulateR()

//DUFBRN()
/** Calculate duff fire properties.
 *
 * @par Original Burnup Description:
!c Duff burning rate (ergo, intensity) and duration

 * @param wdf	Duff loading (kg/m^2, aka W sub d)
 * @param dfm	Ratio of moisture mass to dry organic mass / duff fractional moisture (aka R sub M)
 * @param dfi	Duff fire intensity (aka I sub d)		Units?????
 * @param tdf	Burning duration (aka t sub d)			Units?????
 *
 * @returns On return dfi and tdf are returned in parameters.
 *

! History: Modernized original Burnup subroutine.
 */
void DUFBRN(const double wdf, const double dfm, double& dfi, double& tdf)
{
	double ff;//Fractional duff reduction depth from Brown et al. 1985, (aka F in report equation 4)

	dfi = 0.0;
	tdf = 0.0;
	if ((wdf <= 0.0) || (dfm > 1.96))//This limiting duff moisture should be reviewed!!!!!
	{
		return;
	}
	dfi = 11.25 - 4.05 * dfm;
	ff = 0.837 - 0.426 * dfm;
	tdf = 1.e+04 * ff * wdf / (7.5 - 2.7 * dfm);
}

//DufBrnR()

//GETDAT()

//...

//RETRY()

//ARRAYS
/**
 *
 * @par Original Burnup Description:
!c Orders the fuel description arrays according to the paradigm described in
!c subroutine SORTER and computes the interaction matrix xmat from the array
!c elam and the list alone returned from subroutine OVLAPS.
 *
 * The large number of parameters are hard to handle...

 * The following may be updated on return, in/out: [Refine notes!!!!!!]
 * @param wdry		Ovendry mass loading, kg / sq m. [maxno]
 * @param ash		Mineral content, fraction dry mass. [maxno]
 * @param dendry		Ovendry mass density, kg / cu m. [maxno]
 * @param fmois		Moisture content, fraction dry mass. [maxno]
 * @param sigma		Surface to volume ratio, 1 / m. [maxno]
 * @param htval		Low heat of combustion, J / kg. [maxno]
 * @param cheat		Specific heat capacity, (J / K) / kg dry mass. [maxno]
 * @param condry		Thermal conductivity, W / m K, ovendry. [maxno]
 * @param tpig		Ignition temperature, K. [maxno]
 * @param tchar		Char temperature, K. [maxno]
 * Out:
 * @param diam		Initial diameter, m [by interaction pairs]. [maxkl]
 * @param key 		Ordered index list. [maxno]
 * @param work		Workspace array. [maxno]
 * In:
 * @param  ak		Area influence factor [ak parameter].
 * Out:
 * @param elam		Interaction matrix from OVLAPS. [maxno, maxno]
 * @param alone		Noninteraction fraction list from OVLAPS. [maxno]
 * @param xmat		Consolidated interaction matrix. [maxkl]
 * @param wo		Initial dry loading by interaction pairs. [maxkl]
 * In/out:
 * @param parts		Fuel component names / labels. [maxno]
 * Out?
 * @param list		Intermediary for reordering parts name array. [maxno]
 *            		This is passed in but is not initialized prior
 *            		to that.  It doesn't appear that it is used
 *            		after it is returned.  It appears to only be
 *            		used internal to this routine.
 *            		The same seems to be true for elam and alone?
 * In/out:
 * @param  area		Fraction of site area expected to be covered at
 *             							least once by initial planform area of ea size. [maxno]
 * @param  number	The actual number of fuel classes.  If omitted
 *               	this will be determined from the other inputs.
 *
 * @returns On return many arguments may be sorted, updated, or returned. [MORE!!!!!]

! History: Modernized original Burnup subroutine.
! Several arguments have been removed that were present in the original routine.  The number
! argument has been moved and is now optional and is only used in the interactive context.
*/
void ARRAYS(std::vector<double>& wdry, std::vector<double>& ash, std::vector<double>& dendry,
            std::vector<double>& fmois, std::vector<double>& sigma, std::vector<double>& htval,
            std::vector<double>& cheat, std::vector<double>& condry, std::vector<double>& tpig,
            std::vector<double>& tchar, std::vector<double>& diam, std::vector<int>& key,
            std::vector<double>& work, const double ak, std::vector<std::vector<double>>& elam,
            std::vector<double>& alone, std::vector<double>& xmat, std::vector<double>& wo,
            std::vector<std::string>& parts, std::vector<std::string>& list,
            std::vector<double>& area, const int number)
{
	int numFuelTypes;//The actual number of fuel types, explicit or implied.
	//int j, k, kl, kj;Counters
	int j, k;//Counters
	//double diak, wtk;

	//Testing was done to confrim that explicit initialization of locals was not needed here. [in Fotran]

	//Determine the actual number of fuel types:
	if (number > 0)
	{
		numFuelTypes = number;
	}
	else
	{
		numFuelTypes = wdry.size();
	}

	SORTER(sigma, fmois, dendry, key, number);

	//do j = 1, numFuelTypes
	for (int j = 0; j < numFuelTypes; j++)
	{
		k = key[j];
		list[j] = parts[k];
	}
	for (int j = 0; j < numFuelTypes; j++)
	{
		parts[j] = list[j];
	}

	for (int j = 0; j < numFuelTypes; j++)
	{
		k = key[j];
		work[j] = wdry[k];
	}
	for (int j = 0; j < numFuelTypes; j++)
	{
		wdry[j] = work[j];
	}

	for (int j = 0; j < numFuelTypes; j++)
	{
		k = key[j];
		work[j] = ash[k];
	}
	for (int j = 0; j < numFuelTypes; j++)
	{
		ash[j] = work[j];
	}

	for (int j = 0; j < numFuelTypes; j++)
	{
		k = key[j];
		work[j] = htval[k];
	}
	for (int j = 0; j < numFuelTypes; j++)
	{
		htval[j] = work[j];
	}

	for (int j = 0; j < numFuelTypes; j++)
	{
		k = key[j];
		work[j] = cheat[k];
	}
	for (int j = 0; j < numFuelTypes; j++)
	{
		cheat[j] = work[j];
	}

	for (int j = 0; j < numFuelTypes; j++)
	{
		k = key[j];
		work[j] = condry[k];
	}
	for (int j = 0; j < numFuelTypes; j++)
	{
		condry[j] = work[j];
	}

	for (int j = 0; j < numFuelTypes; j++)
	{
		k = key[j];
		work[j] = tpig[k];
	}
	for (int j = 0; j < numFuelTypes; j++)
	{
		tpig[j] =  work[j];
	}

	for (int j = 0; j < numFuelTypes; j++)
	{
		k = key[j];
		work[j] = tchar[k];
	}
	for (int j = 0; j < numFuelTypes; j++)
	{
		tchar[j] = work[j];
	}

	OVLAPS(wdry, sigma, dendry, ak, fmois, xmat, elam, alone, area, number);

	//do k = 1, numFuelTypes
	for (int k = 0; k < numFuelTypes; k++)
	{
		double diak = 4.0 / sigma[k];
		double wtk = wdry[k];

		//Populate the alone/no companion indexes of the arrays:
		int kl = Loc(k, 0);
		diam[kl] = diak;
		xmat[kl] = alone[k];
		wo[kl] = wtk * xmat[kl];

		//Populate the interacting indexes of the arrays:
		//do j = 1, k
		for (int j = 0; j < numFuelTypes; j++)
		{
			int kj = Loc(k, j);
			diam[kj] = diak;
			xmat[kj] = elam[k][j];
			wo[kj] = wtk * xmat[kj];
		}
	}
}

//SORTER
/**
 *
 * @par Original Burnup Description:
!c Sorts fuel element list in order of increasing size (decreasing sigma)
!c For elements with same size order determined on increasing moisture
!c content (fmois). If items have same size and moisture content, order
!c on the basis of increasing mass density (dryden). "number" elements are
!c included in the list, which has a maximum length of "maxno". The integer
!c list: key(j), j = 1, number holds the indices in order, so other
!c fuel parameters can be ordered and associated as necessary.
 *
 * @param sigma(:)		Surface to volume ratio, 1 / m. [maxno]
 * @param fmois(:)		Moisture content, fraction dry mass. [maxno]
 * @param dryden(:)		Ovendry mass density, kg / cu m. [maxno]
 * @param key(:) 		Ordered index list. [maxno]					//out rather than inout?????
 * @param number		The actual number of fuel classes.  If omitted
 *              		this will be determined from the other inputs.
 * Since SORTER() is only downstream of ARRAYS() making number optional doesn't gain us much.
 *
 * @returns On return sigma, fmois, dryden are reordered and key is set.

! History: Modernized original Burnup subroutine.
! The maxno argument has been removed.  The number argument has been moved and is now optional.
*/
void SORTER(std::vector<double>& sigma, std::vector<double>& fmois, std::vector<double>& dryden,
            std::vector<int>& key, const int number)
{
	int maxNumFuelTypes;		//The maximum number of fuel classes allowed. The input arrays
								//may exceed the actual number. (Not present in original code.)
	int numFuelTypes;			//The actual number of fuel types, explicit or implied.
	//int j, i;					//Counters.
	int i;						//Counter.
	double s, fm, de, keep, usi;			//Hold the values of the current search index. !!!!!
	bool diam, mois, dens, tied;
	bool newIndexFound;//Note: Not present in original code.

	maxNumFuelTypes = sigma.size();//Determine maximum number of fuels from input arrays.
	
	//Determine the actual number of fuel types:
	//if (present(number)) then
	if (number > 0)
	{
		numFuelTypes = number;
	}
	else
	{
		numFuelTypes = maxNumFuelTypes;
	}

	newIndexFound = false;

	//do j = 1, maxNumFuelTypes
	for (int j = 0; j < maxNumFuelTypes; j++)
	{
		key[j] = j;
	}

	//!c Replacement sort: order on increasing size, moisture, density
	//do j = 2, numFuelTypes
	for (int j = 1; j < numFuelTypes; j++)
	{
		//Store the values for this fuel index:
		s = 1.0 / sigma[j];
		fm = fmois[j];
		de = dryden[j];
		keep = key[j];

		//Compare this index (j) with every index before it:
		//do i = (j - 1), 1, -1
		//for (int i = (j - 1); i >= 0; i--)
		for (i = (j - 1); i >= 0; i--)
		{
			usi = 1.0 / sigma[i];
			diam = (usi < s);

			if (diam)
			{
				newIndexFound = true;
				//exit
				break;
			}

			tied = (usi == s);

			if (tied)
			{
				mois = (fmois[i] < fm);
				if (mois)
				{
					newIndexFound = true;
					//exit
					break;
				}

				tied = (fmois[i] == fm);
				if (tied)
				{
					dens = (dryden [i] < de);
					if (dens)
					{
						newIndexFound = true;
						//exit
						break;
					}
				}
			}

			//i is greater than j.
			//Move entry i (the entry we are comparing to) down one:
			sigma[i + 1] = sigma[i];
			fmois[i + 1] = fmois[i];
			dryden[i + 1] = dryden[i];
			key[i + 1] = key[i];
		}

		if (newIndexFound)
		{
			//If a new location has been identified record entry j at i + 1 (below).
			newIndexFound = false;//Reinitialize for the next comparison.
		}
		else
		{
			//i = 0;//Reseting i to 0 will move this entry (j) to postion 1.
			i = -1;//Reseting i to -1 will move this entry (j) to postion 0.
		}

		//Record the values for fuel index j at the identified index:
		sigma[i + 1] = 1.0 / s;
		fmois[i + 1] = fm;
		dryden[i + 1] = de;
		key[i + 1] = keep;
	}
}

//OVLAPS
/**
 *
 * @par Original Burnup Description:
!c Computes the interaction matrix elam(j, k) which apportions the
!c influence of smaller and equal size pieces on eaqh size class for the
!c purpose of establishing the rates at which the elemnts burn out.
!c Input quantities are dryld, the ovendry mass per unit area of each
!c element available for burning after the passage of the igniting surface
!c fire; sigma, the particle's surface/ volume ratio, and dryden, the
!c ovendry mass density of the particle; ak a dimensionless parameter that
!c scales the planform area of a particle to its area of influence. There
!c are "number" separate particle classes, of a maximum number = maxno.
!c It is assumed that the 1ist are ordered on size class (nonincreasing
!c surface/ volume ratio). List "alone" gives the fraction of each loading
!c that is not influenced by any other category.
!
! History: Modernized original Burnup subroutine.
! We modify the original behavior such that a negative value for ak indicates the the value of
! ak / K_a should be calculated.  This requires fmois to be passed in, which was not one of the
! original arguments.
! Several arguments have been removed that were present in the original routine.  The number
! argument has been moved and is now optional.

 * @param dryld		Ovendry mass per unit area of each element (kg/sq m) (= wdry, ...). [maxno]
 * @param sigma		Surface to volume ratio, 1 / m. [maxno]
 * @param dryden	Ovendry mass density, kg / cu m (elsewhere dendry). [maxno]
 * @param ak		Area influence factor (ak / K_a parameter).
 * @param fmois		Moisture fraction of component. [maxno]
 * @param beta		Consolidated interaction matrix (returned). (elsewhere = xmat). [maxkl]
 * @param elam(:,:)	Interaction matrix (returned). [maxno, maxno]
 * @param alone		Non-interacting fraction for each fuel class (returned). [maxno]
 * @param area		Fraction of site area expected to be covered at
 *            		least once by initial planform area of ea size (returned). [maxno]
 * @param number	The actual number of fuel classes.  If omitted
 *              	this will be determined from the other inputs.
 * Since OVLAPS() is only downstream of ARRAYS() making number optional doesn't gain us much.
 *
 * @returns On return beta, elam, alone, and area are returned in parameters.
 */
void OVLAPS(const std::vector<double> dryld, const std::vector<double> sigma,
            const std::vector<double> dryden, const double ak, const std::vector<double> fmois,
            std::vector<double>& beta, std::vector<std::vector<double>>& elam,
            std::vector<double>& alone, std::vector<double>& area, const int number)
{
	double pi;//Convert to a constant?
	int numFuelTypes;//The actual number of fuel types, explicit or implied.
	int j, k, l, kj, kl;//Counters
	double siga;//K_a * diameter
	double a;
	double bb;
	double frac;//Intermediate for calculating alone().
	double K_a;//Value of K_a (ak) parameter, fixed or calculated depending on the mode.

	pi = abs(acos(-1.0));//Calculate pi.

	//Determine the actual number of fuel types:
	//if (present(number)) then
	if (number > 0)
	{
		numFuelTypes = number;
	}
	else
	{
		numFuelTypes = dryld.size();
	}

	//Initialize arrays to 0:
// 	do j = 1, numFuelTypes
// 		alone(j) = 0.0
// 		do k = 1, j
// 			kj = Loc(j, k)
// 			beta(kj) = 0.0
// 		end do
// 		do k = 1, numFuelTypes
// 			elam(j, k) = 0.0
// 		end do
// 	end do
	std::fill(alone.begin(), alone.end(), 0);
	std::fill(beta.begin(), beta.end(), 0);
	for (auto& row : elam)
	{
		std::fill(row.begin(), row.end(), 0);
	}
	//We should allow the vectors of any size, including empty vectors, to be passed in!!!!!

	//do k = 1, numFuelTypes
	//for (int k = 1; k <= numFuelTypes, k++)
	//The indexing for all inputs and outputs are 0 based but we have to convert for Loc():
	for (int k = 0; k < numFuelTypes; k++)//Use 0 indexing.
	{
		//do l = 1, k
		for (int l = 0; l <= k; l++)//Use 0 indexing.
		{
			if (ak > 0.0)//or .ge. ?
			{
				//If a valid value has been supplied use a fixed K_a as in the original Burnup:
				K_a = ak;
			}
			else
			{
				//A negative value indicates that the K_a should be calculated:
				//Calculate ak from the fuel moisture of the smaller or similar fuel member (l):
				//Albini & Reinhardt 1997 Equation 4: K_a = K exp(-B * M^2)
				K_a = 3.25 * exp(-20 * pow(fmois[l], 2));
			}

			//SAV / pi = diameter (units are carried by ak):
			siga = K_a * sigma[k] / pi;

			//kl = Loc(k, l)
			kl = Loc(k + 1, l + 1);//Convert to 1 based indexing for Loc().
			a = siga * dryld[l] / dryden[l];//siga * ? units in meters
			if (k == l)
			{
				bb = 1.0 - exp(-a);			// JMR: FOFEM suggests this can hit 0?
				area[k] = bb;
			}
			else
			{
				bb = std::min(1.0, a);
			}
			beta[kl] = bb;
		}
	}

	// If there is only one fuel type:
	if (numFuelTypes == 1)
	{
		//elam(1, 1) = beta(2)
		elam[0][0] = beta[1];
		alone[0] = 1.0 - elam[0][0];
		return;
		//break;
	}

	//do k = 1, numFuelTypes
	for (int k = 0; k < numFuelTypes; k++)//Use 0 indexing.
	{
		//These inner loops could be combined to simplify the logic and make it more readable!!!!!
		
		frac = 0.0;
		//do l = 1, k
		for (int l = 0; l <= k; l++)//Use 0 indexing.
		{
			//kl = Loc(k, l);
			kl = Loc(k + 1, l + 1);//Convert to 1 based indexing for Loc().
			frac = frac + beta[kl];
		}

		if (frac > 1.0)
		{
			//do l = 1, k
			for (int l = 0; l <= k; l++)//Use 0 indexing.
			{
				//kl = Loc(k, l)
				kl = Loc(k + 1, l + 1);//Convert to 1 based indexing for Loc().
				elam[k][l] = beta[kl] / frac;
			}
			alone[k] = 0.0;
		}
		else
		{
			for (int l = 0; l <= k; l++)//Use 0 indexing.
			{
				//kl = Loc(k, l)
				kl = Loc(k + 1, l + 1);//Convert to 1 based indexing for Loc().
				elam[k][l] = beta[kl];
			}
			alone[k] = 1.0 - frac;
		}
	}
}

//START
/**
 *
 * @par Original Burnup Description:
!c This routine initializes variables prior to starting sequence of calls
!c to subroutine STEP.  On input here, fi is area intensity of spreading
!c fire, dt is the residence time for the spreading fire.  Only physical
!c parameters specified are the fuel array descriptors. To call STEP,
!c one must initialize the following variables.

! Arguments: (by category, not argument order)
! JMR_NOTE: The original comments imply that alfa, diam, and wo should all be intent(in).
! However the code is not consistant with that.
 * @param dt		Spreading fire residence time (s) (= ti, tis, or time elsewhere).
 * @param now		Index marks end of time step.
 * @param wo		Current ovendry loading for the larger of each component pair, kg / sq m.
 *          		Updated on return. [maxkl]
 * @param alfa		Dry thermal diffusivity of component, sq m / s. [maxno]
 * @param dendry	Ovendry density of component, kg / cu m. [maxno]
 * @param fmois		Moisture fraction of component. [maxno]
 * @param cheat		Specific heat capacity of component, J / kg K. [maxno]
 * @param condry	Ovendry thermal conductivity, W / sq m K. [maxno]
 * @param diam		Current diameter of the larger of each
 *            		fuel component pair, m.  Updated on return. [maxkl]
 * @param tpig		Ignition temperature (K), by component. [maxno]
 * @param tchar		tchar = end - pyrolysis temperature (K), by component. [maxno]
 * @param xmat		Table of influence fractions between components. [maxkl]
 * @param tpamb		Ambient temperature (K).
 * @param fi		Current fire intensity (site avg), kW / sq m.
 *
 * Parameters updated (input and output):
 * @param ncalls	Counter of calls to this routine:
 *              	= 0 on first call or reset,
 *              	cumulates after first call.
 *              	JMR_NOTE: This is a strange argument as it is only
 *              	initialized to zero here and is not used.  It is
 *              	returned and passed on to STEP().  It is probably
 *              	initialized here to prevent isses related to
 *              	persistance should more than one simulation be run in
 *              	an interactive session.
 * @param flit		Fraction of each component currently alight. [maxno]
 * @param fout		Fraction of each component currently gone out. [maxno]
 * @param tdry		Time of drying start of the larger of each fuel component pair. [maxkl]	[Units = s?????]
 * @param tign		Ignition time for the larger of each fuel component pair. [maxkl]		[Units = s?????]
 * @param tout		Burnout time of larger component of pairs. [maxkl]
 * @param qcum		Cumulative heat input to larger of pair, J / sq m. [maxkl]
 * @param tcum		Cumulative temp integral for qcum (drying). [maxkl]
 * @param acum		Heat pulse area for historical rate averaging. [maxkl]
 * @param qdot		History (post ignite) of heat transfer rate
 *            		to the larger of each component pair, W / sq m. [maxkl, mxstep]
 * @param ddot		Diameter reduction rate, larger of pair, m / s. [maxkl]
 * @param wodot		Dry loading loss rate for larger of pair. [maxkl]
 * @param work		Workspace array. [maxno]		[Or alternative description in STEP()?????]
 *
 * Constant parameters:
 * @param u			Mean horizontal windspeed at top of fuelbed (m/s).
 * @param d			Fuelbed depth (m).
 * @param r0		Minimum value of mixing parameter.
 * @param dr		Max - min value of mixing parameter.
 *
 * Interactive context:
 * @param number	The actual number of fuel classes.  If omitted
 *              	this will be determined from the other inputs.
 *              	Here the presence of number is used to determine
 *              	if we are in an interactive session as well.

! The constants ch2o and tpdry were included as arguments in the original code.  They have
! chnaged to globals.
! The original comments include hvap as a constant, but is not actually used:
! hvap = heat of vaporization of water J / kg

! History: Modernized original Burnup subroutine.
! Several arguments have been removed that were present in the original routine.  The number
! argument has been moved and is now optional and is only used in the interactive context.
 */
void START(const double dt, const int now, std::vector<double>& wo, std::vector<double>& alfa,
           const std::vector<double> dendry, const std::vector<double> fmois,
           const std::vector<double> cheat, const std::vector<double> condry,
           std::vector<double>& diam, const std::vector<double> tpig,
           const std::vector<double> tchar, const std::vector<double> xmat, const double tpamb,
           const double fi, std::vector<double>& flit, std::vector<double>& fout,
           std::vector<double>& tdry, std::vector<double>& tign, std::vector<double>& tout,
           std::vector<double>& qcum, std::vector<double>& tcum, std::vector<double>& acum,
           std::vector<std::vector<double>>& qdot, std::vector<double>& ddot,
           std::vector<double>& wodot, std::vector<double>& work, const double u, const double d,
           const double r0, const double dr, int& ncalls, const int number)
{
	int numFuelTypes;	//The actual number of fuel types, explicit or implied.
	//integer :: k, l, kl ! Counters.
	int kl;				//Triangular matrix index, 0 based.
	double delm;		//Moisture effect on burning rate (scale factor).
	double heatk;		//Burn rate factor.
	double r;			//Dimensionless mixing parameter.
	double tf;			//Fire environment temperature.
	double ts;			//The charing temperature for a single element (K).
	double thd, tx;
	double dia;			//Diameter for a single element.
	double cpwet;		//Wet specific heat of a single element, J / kg K.
	double fac;			//A factor (radius squared / wet thermal diffusivity) used to convert
	           			//the output of DRYTIM() from dimensionless to actual time.
	double dryt;		//Drying time for a single element.
	double tsd;
	double c;			//Thermal conductivity for a single element.
	double tigk;		//Ignition temperature for a single element.
	double en, e;		//Modified Nusselt number obtained from HEATX() (in different places).
	double trt;			//Minimum ignition time across all fuels (initial estimate?).
	double nlit	;		//Number of elements lit.
	double factor;		//Moisture factor for a single element.
	double hb;			//"Effective" film heat transfer coefficient returned from HEATX().
	double hf;			//Film heat transfer coefficient returned from HEATX().
	double dtign;		//Time to piloted ignition returned from TIGNIT().
	double conwet;		//Wet thermal conductivity for a single element, W / sq m K.
	double aint;
	double ddt;			//Timestep to calculate.  May be less that dt if fuel burns out sooner.
	double dnext;		//Diameter after combustion in this timestep.
	double wnext;		//wo after combustion in this timestep.
	double df;			//

	//Determine the actual number of fuel types:
	//if (present(number)) then
	if (number > 0)
	{
		numFuelTypes = number;
	}
	else
	{
		numFuelTypes = alfa.size();
	}

	/*!c Initialize time varying quantities and, set up work(k)
	!c The diameter reduction rate of fuel component k is given
	!c by the product of the rate of heat tranefer to it, per
	!c unit surface area, and the quantity work(k)*/

	//do k = 1, numFuelTypes
	//for (int k = 0; k < numFuelTypes; k++)//k in 0 based
	for (int k0 = 0; k0 < numFuelTypes; k0++)//k0 = k as base 0 index
	{
		fout[k0] = 0.0;
		flit[k0] = 0.0;
		alfa[k0] = condry[k0] / (dendry[k0] * cheat[k0]);
		//!c effect of moisture content on burning rate (scale factor)
		delm = 1.67 * fmois[k0];
		//!c effect of component mass density (empirical)
		heatk = dendry[k0] / 446.0;
		//!c empirical burn rate factor, J / cu m - K
		heatk = heatk * 2.01e+06 * (1.0 + delm);
		//!c normalize out driving temperature difference (Tfire - Tchar)
		//!c to average value of lab experiments used to find above constants
		work[k0] = 1.0 / (255.0 * heatk);

		//do l = 0, k
		for (int l = 0; l <= k0 + 1; l++)//l in kl space, 0 based
		{
			kl = Loc(k0 + 1, l);
			tout[kl] = rindef;
			tign[kl] = rindef;
			tdry[kl] = rindef;
			tcum[kl] = 0.0;
			qcum[kl] = 0.0;
		}
	}

	//The original code does not initialize acum and qdot.  Failure to do so leads to small
	//variations between runs and changes to results.
	//acum = 0.0;
	std::fill(acum.begin(), acum.end(), 0);
	//qdot = 0.0;
	for (auto& row : qdot)
	{
		std::fill(row.begin(), row.end(), 0);
	}
	//There are a large number of other locals that are not explictly initialized.  Testing was
	//done to confrim that explicit initialization was not needed.

	/*!c Make first estimate of drying start times for all components
	!c These times are usually brief and make little or no difference*/

	r = r0 + 0.25 * dr;
	tf = TEMPF(fi, r, tpamb);
	ts = tpamb;
	if (tf < (tpdry + 10.0))
	{
		//if (present(number)) then ! In an interactive session, preserve the original behavior:
		if (number > 0)//In an interactive session, preserve the original behavior:
		{
			//stop ' Igniting fire cannot dry fuel'
			Stop("Igniting fire cannot dry fuel");
		}
		else//Otherwise signal the condition and return
		{
			//tign = -2.0;//Value signals fuel did not dry.
			std::fill(tign.begin(), tign.end(), -2.0);//Value signals fuel did not dry.
			//write(*, *) "Igniting fire cannot dry fuel."
			Msg.Log("Igniting fire cannot dry fuel");
			return;
		}
	}
	thd = (tpdry - ts) / (tf - ts);
	tx = 0.5 * (ts + tpdry);

	//do k = 1, numFuelTypes
	for (int k0 = 0; k0 < numFuelTypes; k0++)//k0 = k as base 0 index
	{
		factor = dendry[k0] * fmois[k0];
		conwet = condry[k0] + 4.27e-04 * factor;
		//do l = 0, k
		for (int l = 0; l <= k0 + 1; l++)//l in kl space, 0 based
		{
			kl = Loc(k0 + 1, l);
			dia = diam[kl];
			HEATX(u, d, dia, tf, tx, hf, hb, conwet, en);
			//DRYTIM(en, thd, dryt);
			dryt = DRYTIM(en, thd);
			cpwet = cheat[k0] + fmois[k0] * ch2o;
			fac = pow((0.5 * dia), 2) / conwet;
			fac = fac * dendry[k0] * cpwet;
			dryt = fac * dryt;
			tdry[kl] = dryt;
		}
	}

	//!c Next, determine which components are alight in spreading fire

	tsd = tpdry;

	//do k = 1, numFuelTypes
	for (int k0 = 0; k0 < numFuelTypes; k0++)//k0 = k as base 0 index
	{
		c = condry[k0];
		tigk = tpig[k0];
		//do l = 0, k
		for (int l = 0; l <= k0 + 1; l++)//l in kl space, 0 based
		{
			kl = Loc(k0 + 1, l);
			dryt = tdry[kl];
			if (dryt < dt)
			{
				dia = diam[kl];
				ts = 0.5 * (tsd + tigk);
				HEATX(u, d, dia, tf, ts, hf, hb, c, e);
				tcum[kl] = std::max((tf - ts) * (dt - dryt), 0.0);
				qcum[kl] = hb * tcum[kl];
				if (tf > (tigk + 10.0))
				{
					//TIGNIT(tpamb, tpdry, tpig[k0], tf, condry[k0], cheat[k0], fmois[k0], dendry[k0],
					//       hb, dtign);
					dtign = TIGNIT(tpamb, tpdry, tpig[k0], tf, condry[k0], cheat[k0], fmois[k0],
					               dendry[k0], hb);
					trt = dryt + dtign;
					tign[kl] = 0.5 * trt;

					if (dt > trt)
					{
						flit[k0] = flit[k0] + xmat[kl];
					}
				}
			}
		}
	}

	nlit = 0;
	trt = rindef;

	//!c Determine minimum ignition time and verify ignition exists

	//do k = 1, numFuelTypes
	for (int k0 = 0; k0 < numFuelTypes; k0++)//k0 = k as base 0 index
	{
		if (flit[k0] < 0.0)
		{
			nlit = nlit + 1;
		}

		//do l = 0, k
		for (int l = 0; l <= k0 + 1; l++)//l in kl space, 0 based
		{
			kl = Loc(k0 + 1, l);
			trt = std::min(trt, tign[kl]);
		}
	}

	if (nlit == 0)
	{
		//if (present(number)) then ! In an interactive session, preserve the original behavior:
		if (number > 0)//In an interactive session, preserve the original behavior:
		{
			//stop ' START ignites no fuel'
			Stop("START ignites no fuel");
		}
		else//Otherwise signal the condition and return:
		{
			//tign = -1.0;//Value signals fuel did not ignite.
			std::fill(tign.begin(), tign.end(), -1.0);//Value signals fuel did not ignite.
			//write(*, *) "Igniting fire cannot ignite fuel."
			Msg.Log("Igniting fire cannot ignite fuel.");
			return;
		}
	}

	//!c Deduct trt from all time estimates, resetting time origin

	//do k = 1, numFuelTypes
	for (int k0 = 0; k0 < numFuelTypes; k0++)//k0 = k as base 0 index
	{
		//do l = 0, k
		for (int l = 0; l <= k0 + 1; l++)//l in kl space, 0 based
		{
			kl = Loc(k0 + 1, l);
			if (tdry[kl] < rindef)
			{
				tdry[kl] = tdry[kl] - trt;
			}
			if (tign[kl] < rindef)
			{
				tign[kl] = tign[kl] - trt;
			}
		}
	}

	//!c Now go through all component pairs and establish burning rates
	//!c for all the components that are ignited; extrapolate to end time dt

	//do k = 1, numFuelTypes
	for (int k0 = 0; k0 < numFuelTypes; k0++)//k0 = k as base 0 index
	{
		if (flit[k0] == 0.0)
		{
			//do l = 0, k
			for (int l = 0; l <= k0 + 1; l++)//l in kl space, 0 based
			{
				kl = Loc(k0 + 1, l);
				ddot[kl] = 0.0;
				tout[kl] = rindef;
				wodot[kl] = 0.0;
			}
		}
		else
		{
			ts = tchar[k0];
			c = condry[k0];
			//do l = 0, k
			for (int l = 0; l <= k0 + 1; l++)//l in kl space, 0 based
			{
				kl = Loc(k0 + 1, l);
				dia = diam[kl];
				HEATX(u, d, dia, tf, ts, hf, hb, c, e);
				qdot[kl][now] = hb * std::max((tf - ts), 0.0);
				aint = pow((c / hb), 2);
				ddt = dt - tign[kl];
				acum[kl] = aint * ddt;
				ddot[kl] = qdot[kl][now] * work[k0];
				tout[kl] = dia / ddot[kl];
				dnext = std::max(0.0, (dia - ddt * ddot[kl]));
				wnext = wo[kl] * pow((dnext / dia), 2);
				wodot[kl] = (wo[kl] - wnext) / ddt;

				diam[kl] = dnext;
				wo[kl] = wnext;
				df = 0.0;
				if (dnext < 0.0)
				{
					df = xmat[kl];
					wodot[kl] = 0.0;
					ddot[kl] = 0.0;
				}
				flit[k0] = flit[k0] - df;
				fout[k0] = fout[k0] + df;
			}
		}
	}

	ncalls = 0;
}


//FIRINT
/**
 *
 * @par Original Burnup Description:
!c Computes fi = site avg fire intensity given the burning rates of all
!c interacting pairs of fuel components [ wodot ], the mineral ash content
!c of each component [ ash ], the heat of combustion value [ htval ] for
!c each, and the number of fuel components [ number ], where max - maxno.
!c fi is in kW / sq m, while htval is in J / kg.
!
!c fint(k) is the correction to fi to adjust
!c the intensity level to be the local value where size k is burning.

! History: Modernized original Burnup subroutine.
! Several arguments have been removed that were present in the original routine.

 * @param wodot		Burning rates of interacting pairs of fuel components. [maxkl]
 * @param ash		Mineral content, fraction dry mass. [maxno]
 * @param htval		Low heat of combustion, J / kg. [maxno]
 * @param number	The actual number of fuel classes.
 * @param area		Fraction of site area expected to be covered at
 *            		least once by initial planform area of ea size. [maxno]
 * @param fint		A vector to hold the corrected local fire intensity for each fuel type (returned). [maxno]
 * @param fi		Site avg fire intensity (kW / sq m) (returned).
 *
 * @returns fint[] and fi are returned it the parameters.
 */
void FIRINT(const std::vector<double> wodot, const std::vector<double> ash,
            const std::vector<double> htval, const int number, const std::vector<double> area,
            std::vector<double>& fint, double& fi)
{
	double sum;//Running total for fi.
	//integer :: k, l, kl;//Counters
	double wdotk;
	double term;
	double ark;//Area for element k.

	//Constants:
	const double small = 1.0e-06;

	//We could check that the inputs are the same size?

	//Don't assume that fint is the right size!

	sum = 0.0;
	//do k = 1, number
	for (int k = 1; k <= number; k++)//The k fuel types start at 1.
	//for (int k = 0; k < number; k++)
	{
		wdotk = 0.0;
		//do l = 0, k
		for (int l = 0; l <= k; l++)
		{
			int kl = Loc(k, l);
			wdotk = wdotk + wodot[kl];
		}
		term = (1.0 - ash[k]) * htval[k] * wdotk * 1.0e-03;
		ark = area[k];
		if (ark < small)
		{
			fint[k] = term / ark - term;
		}
		else
		{
			fint[k] = 0.0;
		}
		sum = sum + term;
	}

	fi = sum;
}

//STASH

//TIGNIT
/** This routine computes the halfspace surface ignition time under steady radiant heating with
* surface film cooling.
*
! History: Modernized original Burnup subroutine.
!
 *
 * @param tpam	Ambient temperature, K.
 * @param tpdr	Fuel temperature at start of drying, K.
 *            	Currently this is always tpdry, so this argument could be cut.
 * @param tpig	Fuel surface temperature at ignition, K.
 * @param tpfi	Fire enviroriment temperature, K.
 * @param cond	Fuel ovendry thermal conductivity, W / m K.
 * @param chtd	Fuel ovendry specific heat capacity, J / kg K.
 * @param fmof	Fuel moisture content, fraction dry weight.
 * @param dend	Fuel ovendry density, kg / cu m.
 * @param hbar	Effective film heat transfer coefficient [< HEATX] W / sq m K.
	real*4, intent(out) :: tmig	! Predicted time to piloted ignition, s.
 *
 * @returns tmig Predicted time to piloted ignition, s.
 * 
! JMR_NOTE: Since this only has one return value it could be turned into a function.
 * @note The orignal Fortran routine returns tmig as an argument.  !!!!!
 */
double TIGNIT(const double tpam, const double tpdr, const double tpig, const double tpfi,
              const double cond, const double chtd, const double fmof, const double dend,
              const double hbar)//, tmig)
{
	double b03;
	double xlo, xhi, xav;	//Binary search bounds and middle (average).
	double fav;				//Value of ff() for current search value.
	double beta, conw;
	double dtb;				//The temperature increase required to reach the drying temperature.
	double dti;				//The temperature increase required to reach the ignition temperature.
	double ratio, rhoc;
	double tmig;//Return value: Predicted time to piloted ignition, s.

	//Constants:
	const double a03 = -1.3371565;	//ff() parameter 1
	const double a13 = 0.4653628;	//ff() parameter 2
	const double a23 = -0.1282064;	//ff() parameter 3
	const double pinv = 2.125534;
	const double small = 1.e-06;	//Terminate the search when we are this close.
	const double hvap = 2.177e+06;	//Heat of vaporization of water J/kg.
	const double cpm = 4186.0;
	const double conc = 4.27e-04;

	//!c radiant heating equivalent form gives end condition fixes beta

	b03 = a03 * (tpfi - tpig) / (tpfi - tpam);

	//!c find x that solves ff(x) = 0 ; method is binary search

	xlo = 0.0;
	xhi = 1.0;
	// Testing was done to confrim that explicit initialization of other locals was not needed.

	while (true)
	{
		xav = 0.5 * (xlo + xhi);

		//The original code implements this as a statement function:
		//fav = ff(xav)
		/*Original notes:
		!c approximate function of beta to be solved is ff(x) where
		!c x  =  1 / (1 + p * beta)  { Hastings, Approximations for
		!c digital computers } and we employ      pinv  =  1 / p*/
		fav = b03 + xav * (a13 + xav * (a23 + xav));
		/* Or:
		fav = b03 + xav * (a13 + (a23 * xav + xav * xav))
		fav = b03 + (a13 * xav) + (a23 * xav ** 2) + (xav ** 3)*/

		if (abs(fav) < small)
		{
			//exit
			break;
		}
		else if (fav < 0.0)
		{
			xlo = xav;
		}
		else if (fav < 0.0)
		{
			xhi = xav;
		}
	}

	beta = pinv * (1.0 - xav) / xav;
	conw = cond + conc * dend * fmof;
	dtb = tpdr - tpam;
	dti = tpig - tpam;
	ratio = (hvap + cpm * dtb) / (chtd * dti);
	rhoc = dend * chtd * (1.0 + fmof * ratio);
	//tmig = ((beta / hbar) ** 2) * conw * rhoc;
	tmig = (pow((beta / hbar), 2)) * conw * rhoc;

	return tmig;
}

//DRYTIM
/**
!c Given a Nusselt number (enu, actually Biot number = h D / k)
!c and the dimensionless temperature rise required for the start
!c of surface drying (theta), returns the dimerisionless time (tau)
!c needed to achieve it. The time is given multiplied by thermal
!c diffusivity and divided by radius squared. Solution by binary search.
!
! History: Modernized original Burnup subroutine.

 * @param enu Biot number.
 * @param theta Temperature rise required for the start of moisture loss.
			real, intent(out) :: tau	! Time required for the start of moisture loss.

 @note The original Fortran routine returned tau via an argument.
 */
double DRYTIM(const double enu, const double theta)//, tau)
{
	double xl, xh, xm;//The binary search low and high bounds, and the search center.
	double x;
	double approx;
	//int n;//Counter
	double tau;//Return value: Time required for the start of moisture loss.

	//Constants:
	const double p = 0.47047;

	xl = 0.0;
	xh = 1.0;

	/*No rationale is given for the use of 15 cycles in the binary search.  I assume that is was
	determined to be sufficient empirically.  Using a limit for convergence might be better.
	Additionally, since the search middle is calculated at the start ot the loop only 14
	cycles actually inform the result.*/
	//do n = 1, 15
	for (int n = 1; n <= 15; n++)
	{
		xm = 0.5 * (xl + xh);
		approx = ErrorApprox(xm, theta);
		if (approx < 0.0)//or if (ErrorAppox(xm) < 0.0) then
		{
			xl = xm;
		}
		else
		{
			xh = xm;
		}
	}

	x = (1.0 / xm - 1.0) / p;
	//tau = (0.5 * x / enu) **2
	tau = pow((0.5 * x / enu), 2);

	return tau;
}


//HEATX
/** This routine calculates how heat is transfered from the fire environment to a given fuel type.
 *
 * @par Original Burnup Description:
 !c Given horizontal windspeed u at height d [top of fuelbed], cylindrical
!c fuel particle diameter dia, fire environment temperature tf, and mean
!c surface temperature, ts, subroutine returns film heat transfer coefficient
!c hfm and an "effective" film heat transfer coefficient including radiation
!c heat transfer, hbar.  Using the wood's thermal conductivity, cond, the
!c modified Nusselt number [ en ] used to estimate onset of surface drying
!c is returned as well.
!
! History: Modernized original Burnup subroutine.
 *
 * @param u		Mean horizontal windspeed at top of fuelbed (m/s).
 * @param d		Fuelbed depth (m).
 * @param dia	Fuel diameter.									Units?????
 * @param tf	Fire environment temperature.					Units?????
 * @param ts	Mean surface temperature.						Units?????
 * @param hfm	Film heat transfer coefficient (returned).
 * @param hbar	"Effective" film heat transfer coefficient (returned).
 * @param cond	Wood thermal conductivity.						Units?????
 * @param  en	Modified Nusselt number (returned).
 *
 * @returns hfm, hbar, and en are returned via parameters.
 *
 * The arguments are in the original order but it might be good to reorder so the return values are together!!!!!
 */
void HEATX(const double u, const double d, const double dia, const double tf, const double ts,
           double& hfm, double& hbar, const double cond, double& en)
{
	double v;		//Estimate of relative vertical air velocity over fuel element.
	double re;		//Reynolds number (air).
	double enuair;	//Nusselt number.
	double conair;	//(Forced) convection of air?
	double fac;
	double hfmin;	//Film heat transfer coefficient for natural convection (used as minimum value).
	double hrad;	//Radiation contribution.

	//Constants:
	const double g = 9.8;
	const double vis = 7.5e-05;//Kinematic viscosity of hot air.
	const double a = 8.75e-03;
	const double b = 5.75e-05;
	const double rad = 5.67e-08;//Stefan-Boltzmann radiation constant(W/m^2-K^4).
	const double fmfac = 0.382;
	const double hradf = 0.5;//View factor emissivity.

	hfm = 0.0;
	//Testing was done to confrim that explicit initialization of other locals was not needed.

	if (dia > b)
	{
		v = sqrt(u * u + 0.53 * g * d);
		re = v * dia / vis;
		//enuair = 0.344 * (re**0.56)
		enuair = 0.344 * pow(re, 0.56);
		conair = a + b * tf;
		fac = sqrt(abs(tf - ts) / dia);
		hfmin = fmfac * sqrt(fac);
		hfm = std::max((enuair * conair / dia), hfmin);
	}

	hrad = hradf * rad * (tf + ts) * (tf * tf + ts * ts);
	hbar = hfm + hrad;
	en = hbar * dia / cond;
}

//TEMPF
/**
!c Returns a fire environment temperature, TEMPF, given the fire intensity
!c q in kW / square meter, the ambient temperature tamb in Kelvins, and the
!c dimensionless mixing parameter r.
!
! History: Modernized original Burnup subroutine.

 * @param q Fire intensity (kW/m^2).				Units -> Fortran!!!!
 * @param r Dimensionless mixing parameter.
 * @param tamb Ambient temperature (K).
 */
double TEMPF(const double q, const double r, const double tamb)
{
	double term, rlast, den, rnext;
	double tempf;//Return value

	//Constants:
	const double err = 1.0e-04;
	const double aa = 20.0;

	//Testing was done to confrim that explicit initialization of locals was not needed here.

	term = r / (aa * q);
	rlast = r;

	//do
	while (true)
	{
		den = 1.0 + term * (rlast + 1.0) * (rlast * rlast + 1.0);
		rnext = 0.5 * (rlast + 1.0 + r / den);
		if (abs(rnext - rlast) < err)
		{
			tempf = rnext * tamb;
			//return
			break;
		}
		rlast = rnext;
	}

	return tempf;
}

//STEP

/** This routine calculates one timestep of the fuel consumption process.
 *
 * @par Original Burnup Description:
!c Updates status of all fuel component pairs and returns a snapshot


 * Arguments: (by category, not argument order)
 * @param dt		Time step, sec.		[Or: Spreading fire residence time????? Confirm!!!!!]
 * @param now		Index marks end of time step.
 * @param wo		Current ovendry loading for the larger of each component pair, kg / sq m. [maxkl]
 * @param alfa		Dry thermal diffusivity of component, sq m / s. [maxno]
 * @param dendry	Ovendry density of component, kg / cu m. [maxno]
 * @param fmois		Moisture fraction of component. [maxno]
 * @param cheat		Specific heat capacity of component, J / kg K. [maxno]
 * @param condry	Ovendry thermal conductivity, W / sq m K. [maxno]
 * @param diam		Current diameter of the larger of each
 *            		fuel component pair, m.  Updated on return. [maxkl]
 * @param tpig		Ignition temperature (K), by component. [maxno]
 * @param tchar		tchar = end - pyrolysis temperature (K), by component. [maxno]
 * @param xmat		Table of influence fractions between components. [maxkl]
 * @param tpamb		Ambient temperature (K).
 * @param fi		Current fire intensity (site avg), kW / sq m.
 * @param work		Factor of heat transfer rate hbar * (Tfire - Tebar)
 *            		that yields ddot (k). [maxno]
 *
 * Parameters updated (input and output):
 * @param ncalls	Counter of calls to this routine:
 *              	= 0 on first call or reset,
 *              	cumulates after first call.
 * @param flit		Fraction of each component currently alight. [maxno]
 * @param fout		Fraction of each component currently gone out. [maxno]
 * @param tdry		Time of drying start of the larger of each fuel component pair. [maxkl]
 * @param tign		Ignition time for the larger of each fuel component pair. [maxkl]
 * @param tout		Burnout time of larger component of pairs. [maxkl]
 * @param qcum		Cumulative heat input to larger of pair, J / sq m. [maxkl]
 * @param tcum		Cumulative temp integral for qcum (drying). [maxkl]
 * @param acum		Heat pulse area for historical rate averaging. [maxkl]
 * @param qdot		History (post ignite) of heat transfer rate
 *            		to the larger of each component pair, W / sq m. [maxkl, mxstep]
 * @param ddot		Diameter reduction rate, larger of pair, m / s. [maxkl]
 * @param wodot		Dry loading loss rate for larger of pair. [maxkl]
 *
 * Constant parameters:
 * @param u			Mean horizontal windspeed at top of fuelbed (m/s).
 * @param d			Fuelbed depth (m).
 * @param r0		Minimum value of mixing parameter.
 * @param dr		Max - min value of mixing parameter.
 *
 * @param tin		Start of current time step.
 * @param fint		Correction to fi to compute local intensity
 *            		that may be different due to k burning. [maxno]
 * @param fid		Fire intensity due to duff burning ... this is
 *           		used to up the fire intensity for fuel pieces
 *           		that are burning without interacting with others.
 *
 * Interactive context:
 * @param number	The actual number of fuel classes.  If omitted
 *              	this will be determined from the other inputs.

* Differences from the interface of START:
- alfa is input rather than output
- work is nput rather than output
- tin, fint, fid are added.

! History: Modernized original Burnup subroutine.
! Several arguments have been removed that were present in the original routine.  The number
! argument has been moved and is now optional and is only used in the interactive context.
!
! JMR_NOTE: This routine takes a large number of arguments and the order is a bit confusing
! with input and output parameters mixed in the order.
 */
void STEP(const double dt, const int now, std::vector<double>& wo, const std::vector<double> alfa,
          const std::vector<double> dendry, const std::vector<double> fmois,
          const std::vector<double> cheat, const std::vector<double> condry,
          std::vector<double>& diam, const std::vector<double> tpig,
          const std::vector<double> tchar, const std::vector<double> xmat, const double tpamb,
          const double fi, std::vector<double>& flit, std::vector<double>& fout,
          std::vector<double>& tdry, std::vector<double>& tign, std::vector<double>& tout,
          std::vector<double>& qcum, std::vector<double>& tcum, std::vector<double>& acum,
          std::vector<std::vector<double>>& qdot, std::vector<double>& ddot,
          std::vector<double>& wodot, const std::vector<double> work, const double u, const double d,
          const double r0, const double dr, int& ncalls,
          const double tin, const std::vector<double> fint, const double fid, const int number)
{
	/* The constants ch2o and tpdry were included as arguments in the original code.  They have
	! changed to globals.
	! Note: The original code documents the following variable, but it is not actually used.
	!real*4, intent(in) :: hvap			! heat of vaporization of water, J / kg*/

	// Locals: (not in a consistant order)
	int numFuelTypes;//The actual number of fuel types, explicit or implied.
	bool flag;
	double tnow, tnext;//The time of this and the next timestep.
	double tdun;	//The burnout time for a single pair.
	double tgo;		//Time left to burnout.
	double tifi;	//Time when fire ignition phase ended.
	double next;
	double gi;
	int nspan;
	double tst, aint, qqq;
	double tav1, tav2, tav3, tavg;//Time over which to perform averaging.
	double tbar;
	int index;											//Is index a keyword?????
	double qdsum;	//Sum of heat transfer (W/m^2 * s = J/m^2 ?).
	double qdavg;	//Average heat transfer...
	double deltim, rate, dryt, dqdt;
	double qd;
	double dteff, heff, delt;
	double factor;	//Moisture factor for a single element.
	double dtef;
	double he;		//qcum / tcum
	double tf;		//Fire environment temperature.
	double ts;		//The charing temperature for a single element (K).
	double biot;	//Biot number for a single element.
	double cpwet;	//Wet specific heat of a single element, J / kg K.
	double c;		//Thermal conductivity for a single fuel component.
	double conwet;	//Wet thermal conductivity for a single element, W / sq m K.
	double ddt;		//Timestep to calculate.  May be less that dt if fuel burns out sooner.
	double dia;		//Diameter for a single fuel component (kl).
	double dnext;	//Diameter after combustion in this timestep.
	double wnext;	//wo after combustion in this timestep.
	double dtcum;
	double dtlite;	//Time to ignition returned by TIGNIT().
	double e;
	double fac;		//A factor (radius squared / wet thermal diffusivity) used to convert
					//the output of DRYTIM() from dimensionless to actual time.
	double hb, hf;
	double r;		//Dimensionless mixing parameter.
	double tfe;
	double thd;
	double tlit;	//Ignition time for a single fuel component.
	double tspan;
	double dtemp;

	//integer :: k, l, mu, kl; //Counters (kl is a bit different)
	int kl;			//Triangular matrix index, 0 based.

	//There are a large number of locals in this routine that are not explictly initialized.
	//Testing was done to confrim that explicit initialization was not needed.

	//Determine the actual number of fuel types:
	if (number > 0)
	{
		numFuelTypes = number;
	}
	else
	{
		numFuelTypes = alfa.size();
	}

	ncalls = ncalls + 1;
	tnow = tin;
	tnext = tnow + dt;
	//!c tifi = time when fire ignition phase ended (at now = 1)
	//tifi = tnow - float(now - 1) * dt;
	tifi = tnow - static_cast<double>(now - 1) * dt;
	next = now + 1;

	//kLoop : do k = 1, numFuelTypes				//!!!!!!
	for (int k0 = 0; k0 < numFuelTypes; k0++)//k0 = k as base 0 index
	{
		c = condry[k0];
		//lLoop : do l = 0, k
		for (int l = 0; l <= k0; l++)//Start major lLoop.
		{
			//kl = Loc(k, l);
			kl = Loc(k0 + 1, l);
			tdun = tout[kl];

			//!c See if k of (k, l) pair burned out

			if (tnow >= tdun)
			{
				ddot[kl] = 0.0;
				wodot[kl] = 0.0;
				//cycle lLoop				//!!!!!!
				continue;//Back to start of lLoop
			}
			if (tnext >= tdun)
			{
				tgo = tdun - tnow;
				ddot[kl] = diam[kl] / tgo;

				wodot[kl] = wo[kl] / tgo;
				wo[kl] = 0.0;
				diam[kl] = 0.0;
				//cycle lLoop				//!!!!!!
				continue;//Back to start of lLoop
			}

			//!c k has not yet burned out ... see if k of (k, l) pair is ignited

			tlit = tign[kl];
			if (tnow >= tlit)
			{
				ts = tchar[k0];
				if (l == 0)
				{
					r = r0 + 0.5 * dr;
					gi = fi + fid;
				}
				//if ((l != 0) && (l != k))
				if ((l != 0) && (l != k0 + 1))
				{
					r = r0 + 0.5 * (1.0 + flit[l]) * dr;
					gi = fi + fint[k0] + flit[l] * fint[l];
				}								//JMR: Link with following if!!!!!!
				//if (l == k)
				if (l == k0 + 1)
				{
					r = r0 + 0.5 * (1.0 + flit[k0]) * dr;
					gi = fi + flit [k0] * fint[k0];
				}
				tf = TEMPF(gi, r, tpamb);
				dia = diam[kl];
				HEATX(u, d, dia, tf, ts, hf, hb, c, e);
				qqq = hb * std::max((tf - ts), 0.0);
				tst = std::max(tlit, tifi);
				//nspan = std::max(l, nint((tnext - tst) / dt));
				//nspan = std::max(l, std::round((tnext - tst) / dt));
				nspan = std::max(l, static_cast<int>(std::round((tnext - tst) / dt)));//Ugly!!!!!
				if (nspan <= mxstep)
				{
					qdot[kl][nspan] = qqq;
				}
				else//if (nspan > mxstep)
				{
					//do mu = 2, mxstep
					for (int mu = 1; mu < mxstep; mu++)//Is this right?????
					{
						qdot[kl][mu - 1] = qdot[kl][mu];
					}
					qdot[kl][mxstep] = qqq;
				}
				aint = pow((c / hb), 2);
				acum[kl] = acum[kl] + aint * dt;

				//Time over which to perform averaging:
				tav1 = tnext - tlit;//Time since ignition.
				tav2 = acum[kl] / alfa[k0];//Measure of square of distance heat has penetrated fuel.
				tav3 = pow((dia / 4.0), 2) / alfa [k0];//Measure of time heat takes to reach center of fuel.
				//tavg = std::min(tav1, tav2, tav3);
				tavg = std::min({tav1, tav2, tav3});

				index = 1 + std::min(nspan, mxstep);
				qdsum = 0.0;
				tspan = 0.0;
				deltim = dt;

				//Calculate qdsum (sum of heat transfer (W/m^2 * s = J/m^2)):
				//do
				while (true)
				{
					index = index - 1;
					if (index == 1)
					{
						deltim = tnext - tspan - tlit;
					}

					if ((tspan + deltim) >= tavg)
					{
						deltim = tavg - tspan;
					}

					qdsum = qdsum + qdot[kl][index] * deltim;
					tspan = tspan + deltim;

					if ((tspan < tavg) && (index > 1))
					{
						//cycle
						continue;
					}
					else
					{
						//exit
						break;
					}
				}

				qdavg = std::max(qdsum / tspan, 0.0);
				ddot[kl] = qdavg * work[k0];
				dnext = std::max(0.0, dia - dt * ddot[kl]);

				wnext = wo[kl] * pow((dnext / dia), 2);
				if ((dnext == 0.0) && (ddot[kl] > 0.0))
				{
					tout[kl] = tnow + dia / ddot[kl];
				}
				if ((dnext > 0.0) && (dnext < dia))
				{
					rate = dia / (dia - dnext);
					tout[kl] = tnow + rate * dt;
				}
				if (qdavg <= 20.0)
				{
					tout[kl] = 0.5 * (tnow + tnext);
				}
				ddt = std::min(dt, (tout[kl] - tnow));
				wodot[kl] = (wo[kl] - wnext) / ddt;
				diam[kl] = dnext;
				wo[kl] = wnext;
				//cycle lLoop						//!!!!!!
				continue;//Back to start of lLoop
			}//if (tnow >= tlit)

			//!c See if k of (k, l) has reached outer surface drying stage yet

			dryt = tdry[kl];
			if ((tnow >= dryt) && (tnow < tlit))
			{
				if (l == 0)
				{
					r = r0;
					gi = fi + fid;
				}
				//if (l == k)
				if (l == k0 + 1)
				{
					r = r0;
					gi = fi;
				}
				//if ((l != 0) && (l != k))
				if ((l != 0) && (l != k0 + 1))
				{
					r = r0 + 0.5 * flit[l] * dr;
					gi = fi + flit[l] * fint[l];
				}
				tf = TEMPF(gi, r, tpamb);
				ts = tpamb;
				dia = diam[kl];
				HEATX(u, d, dia, tf, ts, hf, hb, c, e);
				dtemp = std::max(0.0, (tf - ts));
				dqdt = hb * dtemp;
				qcum[kl] = qcum[kl] + dqdt * dt;
				tcum[kl] = tcum[kl] + dtemp * dt;
				dteff = tcum[kl] / (tnext - dryt);
				heff = qcum[kl] / tcum[kl];
				tfe = ts + dteff;
				dtlite = rindef;

				if (!(tfe <= (tpig[k0] + 10.0)))
				{
					//TIGNIT(tpamb, tpdry, tpig[k0], tfe, condry[k0], cheat[k0], fmois[k0], dendry[k0],
					//       heff, dtlite);
					dtlite = TIGNIT(tpamb, tpdry, tpig[k0], tfe, condry[k0], cheat[k0], fmois[k0],
					                dendry[k0], heff);
				}
				tign[kl] = 0.5 * (dryt + dtlite);

				//!c If k will ignite before time step over, must interpolate

				if (tnext > tign[kl])
				{
					ts = tchar[k0];
					HEATX(u, d, dia, tf, ts, hf, hb, c, e);

					//qdot(kl, 1) = hb * max((tf - ts), 0.0)
					qdot[kl][0] = hb * std::max((tf - ts), 0.0);
					//qd = qdot(kl, 1)
					qd = qdot[kl][0];
					ddot[kl] = qd * work[k0];
					delt = tnext - tign[kl];
					dnext = std::max(0.0, dia - delt * ddot[kl]);
					wnext = wo[kl] * pow((dnext / dia), 2);
					if (dnext == 0.0) tout[kl] = tnow + dia / ddot[kl];//Oneliner!!!!!
					if ((dnext > 0.0) && (dnext < dia))
					{
						rate = dia / (dia - dnext);
						tout[kl] = tnow + rate * dt;
					}
					if (tout[kl] > now)
					{
						ddt = std::min(dt, (tout[kl] - tnow));
						wodot[kl] = (wo[kl] - wnext) / ddt;
					}
					else
					{
						wodot[kl] = 0.0;
					}
					diam[kl] = dnext;
					wo[kl] = wnext;
				}
				//cycle lLoop					//!!!!!
				continue;//Back to start of lLoop
			}

			//!c If k of (k, l) still coming up to drying temperature, accumulate
			//!c heat input and driving temperature difference, predict drying start

			if (tnow < dryt)
			{
				factor = fmois[k0] * dendry[k0];
				conwet = condry[k0] + 4.27e-04 * factor;
				if (l == 0)
				{
					r = r0;
					gi = fi + fid;
				}
				//if (l == k)
				if (l == k0 + 1)
				{
					r = r0;
					gi = fi;
				}
				//if ((l != 0) && (l != k))
				if ((l != 0) && (l != k0 + 1))
				{
					r = r0 + 0.5 * flit[l] * dr;
					gi = fi + flit[l] * fint[l];
				}
				tf = TEMPF(gi, r, tpamb);
				if (tf <= (tpdry + 10.0))
				{
					//cycle lLoop				//!!!!!!!
					continue;//Back to start of lLoop
				}
				dia = diam[kl];
				ts = 0.5 * (tpamb + tpdry);
				HEATX(u, d, dia, tf, ts, hf, hb, c, e);
				dtcum = std::max((tf - ts) * dt, 0.0);
				tcum[kl] = tcum[kl] + dtcum;
				qcum[kl] = qcum[kl] + hb * dtcum;
				he = qcum[kl] / tcum[kl];
				dtef = tcum[kl] / tnext;
				thd = (tpdry - tpamb) / dtef;
				if (thd > 0.9)
				{
					//cycle lLoop			//!!!!!!
					continue;//Back to start of lLoop
				}
				biot = he * dia / conwet;
				//call DRYTIM(biot, thd, dryt)
				dryt = DRYTIM(biot, thd);

				cpwet = cheat[k0] + ch2o * fmois[k0];
				fac = pow((0.5 * dia), 2) / conwet;
				fac = fac * cpwet * dendry[k0];
				tdry[kl] = fac * dryt;

				if (tdry[kl] < tnext)
				{
					ts = tpdry;
					HEATX(u, d, dia, tf, ts, hf, hb, c, e);
					dqdt = hb * (tf - ts);
					delt = tnext - tdry[kl];
					qcum[kl] = dqdt * delt;
					tcum[kl] = (tf - ts) * delt;
					tbar = 0.5 * (tpdry + tpig[k0]);

					//!c See if ignition to occur before time step complete

					if (tf <= (tpig[k0] + 10.0))
					{ 
						//cycle lLoop					//!!!!!
						continue;//Back to start of lLoop
					}
					//TIGNIT(tpamb, tpdry, tpig[k0], tf, condry[k0], cheat[k0], fmois[k0], dendry[k0],
					//       hb, dtlite);
					dtlite = TIGNIT(tpamb, tpdry, tpig[k0], tf, condry[k0], cheat[k0], fmois[k0],
					                dendry[k0], hb);
					tign[kl] = 0.5 * (tdry[kl] + dtlite);

					if (tnext > tign[kl])
					{
						ts = tchar[k0];
						//qdot(kl, 1) = hb * max((tf - ts), 0.0)
						qdot[kl][1] = hb * std::max((tf - ts), 0.0);
					}
				}
			}
		}//lLoop
	}//kLoop

	//!c Update fractions ignited and burned out, to apply at next step start

	//do k = 1, numFuelTypes
	for (int k0 = 0; k0 < numFuelTypes; k0++)//k0 = k as base 0 index
	{
		flit[k0] = 0.0;
		fout[k0] = 0.0;
		//do l = 0, k
		for (int l = 0; l <= k0; l++)
		{
			kl = Loc(k0 + 1, l);
			flag = (tnext >= tign[kl]);
			if (flag && (tnext <= tout[kl]))
			{
				flit[k0] = flit[k0] + xmat[kl];
			}
			if (tnext > tout[kl])
			{
				fout[k0] = fout[k0] + xmat[kl];
			}
		}
	}
}

//SUMMARY

//Loc
/** This function converts the indexes of the pairwise fuel interaction triangular matrix space
 * to indexes of the arrays used to represent it for several computed variables.
 *
 * @param k Triangular matrix column (row) index, (1 - number of fuel types).
 * @param l Triangular matrix row (column) index, (0 - k), = partner.
 *          This index starts at 0, which represent the "no companion" pairs.
 *
 * @returns The compact array index representing triangular matrix position [k, l]. 
 * 
 * @par Indexes:
 * The matrix indices represent the fuel types and are 1 based, i.e. 1 = fuel type 1.  The l index
 * represent the partner fuel with 0 representing no partner.  This indexing approach makes
 * reasonable sense.  In the orignal Fortran it also alignes with the default 1 based array
 * indexing.  There is no easy way to make this more C-ish since one dimension is already 0 based.
 * Therefore we keep the indexing the same.  Calling code must keep the index model in mind.
 * However, this implementation differes from the Fortran code in that we assume the arrays used to
 * represent triangular matrix data use native 0 indexing.
 * 
 * @note: This will only return valid (occupied) coordinates of the triangular matrix.
 * The code assumes maxno and maxkl are the maximum dimensions.
 * Further error checking would require that the number of fuel classes be know. ?????
 
! History: This function was originally implemented as a statement function defined in seven places
! in the original code.
 
 */
int Loc(const int k, const int l)
{
	int loc;//Return value: Index in a compact array representing the triangular matrix values.

	//Validity checking:
	if ((k < 1) || (k > maxno))
	{
		//print *, "Loc(): Invalid value of k ", k
		//std::cout << "Loc(): Invalid value of k " << k << std::endl;
		Stop("Loc(): Invalid value of k " + std::to_string(k));
	}

	if ((l < 0) || (l > k))
	{
		//print *, "Loc(): Invalid value of l ", l
		//std::cout << "Loc(): Invalid value of l " << l << std::endl;
		Stop("Loc(): Invalid value of l " + std::to_string(l));
	}

	loc = k * (k + 1) / 2 + l;
	loc -= 1;//Convert to 0 based index.

	//if ((loc < 1) || (loc > maxkl))
	if ((loc < 0) || (loc > (maxkl - 1)))
	{
		//print *, "Loc(): Invalid index returned ", loc
		std::cout << "Loc(): Invalid index returned " << loc << std::endl;
		Stop("Loc(): Invalid index returned " + std::to_string(loc));
	}

	return loc;
}

//ErrorApprox
/**
I believe this function is the approximation of the complementary error function as
described in Albini 1995 and Albini & Reinhardt 1995.  It was obtained from Hastings (et
al.) 1955, but I cann't identify the equation in that reference.

History: This function was originally implemented as a statement function in DRYTIM(), as f().

* @param h
& @param theta Temperature rise required for the start of moisture loss.

 */
double ErrorApprox(const double h, const double theta)
{
	double approx;//Return value.

	//Constants:
	const double a = 0.7478556;
	const double b = 0.4653628;
	const double c = 0.1282064;

	approx = h * (b - h * (c - h)) - (1.0 - theta) / a;

	return approx;
}

//AskForReal1

//...

//SaveStateToFile
/** Output the state of the simulation at the current timestep to file:
 * Sequential calls to this routine will produce a full history of the simulated fire.
 *
 * Successful output from this routine is treated as non-critical as it doesn't impact the
 * simulation process.  Output failures are reported but are not treated as fatal.
 *
 * The file name is currently fixed.  A history file from a previous run will prevent a new one
 * from being created.  In the future we may automatically number the file if it already exists.
 * Allowing the file name to be specified would only work in some contexts.  For example, we
 * currently have no way to pass in a file name via the R interface.
 *
 * @param ts		Current timestep count.
 * @param time		Current time (s).
 * @param number	Actual number of fuel components.
 * @param parts		Fuel component names / labels. [maxno]
 *
 * All the outputs from START() and STEP():
 * Currently only some of these are passed in to be saved. In the future others may be added.
 * @param wo		Current ovendry loading for the larger of each component pair, kg / sq m. [maxkl]
 * @param diam		Current diameter of the larger of each fuel component pair, m. [maxkl]	!!!!!
	!real*4, intent(in) :: flit(maxno)		Fraction of each component currently alight.
	!real*4, intent(in) :: fout(maxno)		Fraction of each component currently gone out.
 * The following are recomputed at each time point but only the final values would be needed:
	!real*4, intent(in) :: tdry(maxkl)		! Time of drying start of the larger of each fuel component pair.
	!real*4, intent(in) :: tign(maxkl)		! Ignition time for the larger of each fuel component pair. [maxkl]
	!real*4, intent(in) :: tout(maxkl)		! Burnout time of larger component of pairs.
	!real*4, intent(in) :: qcum(maxkl)		! Cumulative heat input to larger of pair, J / sq m.
	!real*4, intent(in) :: tcum(maxkl)		! Cumulative temp integral for qcum (drying).
	!real*4, intent(in) :: acum(maxkl)		! Heat pulse area for historical rate averaging.
 * This includes time so again only the final value would be needed:
	!real*4, intent(in) :: qdot(maxkl, mxstep)	! History (post ignite) of heat transfer rate
												! to the larger of each component pair, W / sq m..
	!real*4, intent(in) :: ddot(maxkl)		! Diameter reduction rate, larger of pair, m / s.
	!real*4, intent(in) :: wodot(maxkl)		! Dry loading loss rate for larger of pair.
	!real*4, intent(in) :: work(maxno)	! Workspace array.
 * @param fi		Current fire intensity (site avg), kW / sq m.
 *
 * @returns Nothing.
 *
 * @par Format for the variable output:
 * The data is written in long format with tab delimited fields:
 * timestep (integer), time (float), variable (string), value (float), and IDs (strings).
 * The IDs are currently only used to identify the fuel pairs.  If that is the only use they
 * should be renamed.
 *
 ! History: Added for module.
 */
void SaveStateToFile(const int ts, const double time, const int number,
                     const std::vector<std::string> parts, const std::vector<double> wo,
                     const std::vector<double> diam, const double fi)//const
{
	//Local constants:
	//const char histFileName[] = "BurnupHistory.txt";//histFile in Fortran!!!!!
	const std::string histFileName("BurnupHistory.txt");//histFile in Fortran!!!!!
	const char delim = '\t';//Delimiter = tab character

	//Locals:
	//integer :: k, l, kl ! Counters.
	//int kl;//Triangular matrix index, 0 based.

	std::string fuelName;//Name of the (larger) fuel type.
	std::string compName;//The name of the companion/partner fuel.

	std::ofstream histFile;

	if (SaveHistory)
	{
		//Create or open the history file:
		if (ts == 0)//In the first timestep create and set up the file:
		{
			//open(hUnit, file = histFile, status = 'NEW', iostat = openStat, iomsg = ioMsg)
			//std;:iostream histFile(histFileName);
			histFile.open(histFileName);
			//if (openStat .ne. 0) then
			if (!histFile.is_open())
			{
				//print *, "Can't create file: ", histFile, ", Error: ", openStat, ", Message: ", ioMsg
				//The error message will likely include the file name so we can omit it.
				//print warnFmt, "Can't create file, Error: ", openStat, ", Message: ", ioMsg
				
				Warning("Can't create file: " + histFileName + ", Error: " + std::strerror(errno));
				SaveHistory = false;//If we can't create the file don't try anything further.
				return;
			}

			//Write a column header for the file:
			//write(hUnit, '(a)') "Timestep	TimeSec	Variable	Value	ID1	ID2"
			histFile << "Timestep	TimeSec	Variable	Value	ID1	ID2" << std::endl;
		}
		else//Reopen the file and append:
		{
			//open(hUnit, file = histFile, position = 'APPEND', status = 'OLD', &
			//	 iostat = openStat, iomsg = ioMsg)
			histFile.open(histFileName);
			//if (openStat .ne. 0) then
			if (!histFile.is_open())
			{
				//print warnFmt, "Can't reopen file, Error: ", openStat, ", Message: ", ioMsg
				
				Warning("Can't reopen file: " + histFileName + ", Error: " + std::strerror(errno));
				SaveHistory = false;//Assume the error will persist so don't try again.
				return;
			}
		}

		//do k = 1, number
		for (int k0 = 0; k0 < number; k0++)
		{
			fuelName = parts[k0];

			//do l = 0, k
			for (int l = 0; l < number; l++)
			{
				int kl = Loc(k0 + 1, l);//Triangular matrix index, 0 based.

				//Get the name of the partner component:
				if (l == 0)
				{
					compName = noCmpStr;//NEED to make a string for this to work!!!!
				}
				else
				{
					compName = parts[l];
				}

				//Fuel loading:
				//write(hUnit, formatDelim, iostat = writeStat, iomsg = ioMsg) &
				//	  ts, time, "w_o", wo(kl), trim(fuelName), trim(compName)
				histFile << ts << delim << time << delim << "w_o" << delim << wo[kl] << delim
				         << fuelName << delim << compName << '\n';

				//Particle diameter:
				//write(hUnit, formatDelim, iostat = writeStat, iomsg = ioMsg) &
				//	  ts, time, "Diameter", diam(kl), trim(fuelName), trim(compName)
				histFile << ts << delim << time << delim << "Diameter" << delim << diam[kl] << delim
				         << fuelName << delim << compName << '\n';
			}
		}

		//Average fire intensity:
		//write(hUnit, formatDelim, iostat = writeStat, iomsg = ioMsg) &
		//	  ts, time, "FireIntensity", fi, "NA", "NA"
		histFile << ts << delim << time << delim << "FireIntensity" << delim << fi << delim
				         << "NA" << delim << "NA" << std::endl;

		histFile.close();//Close the file.
		//Should check for any write errors here!!!!!
	}//(SaveHistory)
}
