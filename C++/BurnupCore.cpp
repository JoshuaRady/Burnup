/***************************************************************************************************
BurnupCore.cpp
Burnup Wildfire Fuel Consumption Model

Original Fortan code by: Frank A. Albini
Edited, modified, and ported to C++ by: Joshua M. Rady
Woodwell Climate Research Center
Started: 12/30/2024
Reference: Proj. 11 Exp. 22

This is an reimplementation of the Burnup wildfire fuel consumption model in C++.

	This module is based on the Burnup model by Frank Albini and collaborators (see references).
The original Burnup source code (Fortran IV/66/77, fixed form) was previously updated to modern
Fortran (2003+) (Burnup.f90) and was reformulated into a form that can be compiled as a linkable
module or shared library (BurnupMod.f90 & BurnupInteractive.f90).  Here we have ported the linkable
functionality of the later to C++.

	The goal was to provide a version of the Burnup model that could be easily embedded in other
code or coupled with other models.  The C++ code produces identical results to the original Fortran
code in testing and has limit dependancies.  The intent is that calling code only need to know
Burnup's input and outputs.

	A R interface function has been included for the main program entry point to allow the code to
interface with R when compiled as a shared library (.so file).

 	The routines that run the UI of the original interactive program have not been ported.  They may
be ported in the future in a separate file.

Formating:
	I have used tabs for code indenting and for alignment of comments.  I have used 4-space
equivalent tabs for layout purposes here.

In-code Documentation:
	Many comments have been added to the code to increase it's readability and to document changes.
Comments from the original code are marked with '!c' or are labeled.  Functions are documented with
Doxygen comments.

	Routines from the original propgram maintain their all caps format while new routines use
camelcase.

Dependancies:
	This code was written with limited dependancies.  This code relies on the C++11(+?) standard
template library (STL).  It also uses a small number of utilities from the Fireweed Wildfire Code
Library [reference to come].

References:
	The original source code was obtained from:

Frank A. Albini.
Program BURNUP, a simulation model of the burning of large woody natural fuels.
Missoula, MT: USDA Forest Service. Unpublished report. Research Grant INT-92754-GR. 1994.
Appendix B.

	The original report is held by the National Forest Service Library in Fort Collins.

Albini, F.A. and Reinhardt, E.D.
Improved calibration of a large fuel burnout model.
International Journal of Wildland Fire, 7(1): 21-28, 1997.

	This paper documents an alternate way to calculate the ak / K_a parameter.

Caveats:
	This code compiles, runs, and produces output identical to the original Burnup code in tests but
the code is under active development and errors could be present that have not been identified.

	There is no licence provided for the original code.  It is though to be open by provenance, but
that may not be correct.
***************************************************************************************************/

#include <algorithm>//For max(), min(), fill().
#include <cerrno>
#include <cmath>//For pow(), sqrt(), abs(), acos()..
#include <cstring>//For strerror().
#include <fstream>
#include <iostream>//Or just <ostream>?
#include <vector>

#include "BurnupCore.h"
#include "FireweedMessaging.h"
#include "FireweedUtils.h"//Just for SameLengths().

/* Program level dimensional constants:
In the original code the parameters maxno and maxkl were passed into all routines that need them.
In our Fortran module we changed them to module level scope.  This allowed us to reduce the number
of arguments to routines.  Without the interactive UI they are not yet needed here.

The original interactive program has a fixed maximum number of fuel components or types = maxno.
This is used to build fixed size data structures.  The actual number of fuel elements for a
simulation may be less than this maximum so some array indexes may be empty during computations.
The original program fixes this maximum arbitrarily at 10 fuel components.*/
//const int maxno = 12;

/*The maximum number of non-zero entries in the triangular matrix of fuel interaction pairs is
calculated from maxno to size arrays that hold data is this form.*/
//Add one to one dimension for the 'no companion' interaction element:
//const int maxkl = maxno * (maxno + 1) / 2 + maxno;

/*We have removed these fixed size assumptions from our programatic interface to the model, though
they remain in the intactive UI Fortran version.  The number of fuel types is passed in at the start
of the simulation and other functions can infer it from inputs.  However, we still need to maintain
a global record of the number of fuel types for Loc() to be able to validate its inputs.  See
NumFuelTypes below. This is actual number of fuel types in use, unlike maxno.  When Simulate() is
called all vectors will always be full.

In the code that follows vector parameter descriptions contain their sizes in square brackets.
Currently those marked [maxno] will always equal the actual number of fuel types and those marked
[maxkl] the triangular matrix length for that number of fuels.  The labeling might be improved.*/

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

//Globals:------------------------------------------------------------------------------------------
bool SaveHistory = false;//Should fire history be output to file?
int NumFuelTypes = 0;

//Code:---------------------------------------------------------------------------------------------

//Simulate():
/** Perform a simulation with prescribed inputs and return fuel consumption properties.
 * The main fire properties are returned as output parameters.  Optionally an additional detailed
 * fire history can be output to file.
 *
 * Igniting fire and environmental data:
 * @param[in,out] fi	Current fire intensity (site avg), kW / sq m
 * The value passed in for fi is the igniting fire intensity.  The variable is later reused and
 * updated by FIRINT().  It is passed on to other routines that use but do not change it.
 * These two uses could be separated.  The value returned is the final intensity, which might
 * be of use.  A history would be more valuable.
 *
 * @param[in] ti		Igniting fire residence time (s).
 * @param[in] u			Mean horizontal windspeed at top of fuelbed (m/s).
 * @param[in] d			Fuelbed depth (m).
 * @param[in] tpamb		Ambient temperature (K).
 * 
 * Internal and control variables:
 * @param[in] ak		Area influence factor (ak / K_a parameter).
 *              		We modify the original behavior such that a negative value indicates that
 *              		the value of ak / K_a should be calculated according to Albini & Reinhardt
 *              		1997.
 * @param[in] r0		Minimum value of mixing parameter.
 * @param[in] dr		Max - min value of mixing parameter.
 * @param[in,out] dt	Time step for integration of burning rates (s).
 *                  	On completion contains the time the fire went out.
 *                  	A value of -1 indicates the fuel did not ignite.  A value
 *                  	of -2 indicates the fuel did not complete drying.  In such
 *                  	cases most of remaining return variables will not be
 *                  	meaningful.
 * 
 * Considering removing these two.  See below:
 * @param[in] wdf		Duff loading (kg/m^2, aka W sub d).
 * @param[in] dfm		Ratio of moisture mass to dry organic mass /
 *               		duff fractional moisture (aka R sub M).
 * @param[in] ntimes	Number of time steps to run.  Move down?
 * @param[in] number	The number of fuel classes.
 *
 * Fuel component property arrays:  The values will not change but they may be reordered.
 * Returning the reordered arrays may be overkill.  The revised order might be sufficient.
 * However, setting these to inout allows the values to be reordered internally by SORTER(),
 * which eliminates the need for parallel local variables.
 * @param[in,out] parts		Fuel component names / labels. [maxno]
 * @param[in,out] wdry		Ovendry mass loading, kg/sq m. [maxno]
 * @param[in,out] ash		Mineral content, fraction dry mass. [maxno]
 * @param[in,out] htval		Low heat of combustion (AKA heat content), J / kg. [maxno]
 * @param[in,out] fmois		Moisture fraction of component. [maxno]
 * @param[in,out] dendry	Ovendry mass density, kg / cu m. [maxno]
 * @param[in,out] sigma		Surface to volume ratio, 1 / m. [maxno]
 * @param[in,out] cheat		Specific heat capacity, (J / K) / kg dry mass. [maxno]
 * @param[in,out] condry	Thermal conductivity, W / m K, ovendry. [maxno]
 * @param[in,out] tpig		Ignition temperature, K. [maxno]
 * @param[in,out] tchar		Char temperature, K. [maxno]
 *
 * Calculated outputs:
 * The following are the main variables output by SUMMARY(): [name], fr, ti, to, wd, di
 *
 * wo should be moved to the front in any case because it is the most valuable output.  This
 * would have the advantage of making the size() shorthand shorter.
 * JMR_Note: No longer in argument order!!!!!
 * @param[out] wo		Current ovendry loading for the larger of each component pair, kg/sq m. [maxkl]
 * @param[out] xmat		Table of influence fractions between components. [maxkl]
 * @param[out] tign		Ignition time for the larger of each fuel component pair, s. [maxkl]
 * @param[out] tout		Burnout time of larger component of pairs, s. [maxkl]
 * @param[out] diam		Current diameter of the larger of each fuel component pair, m. [maxkl]
 *
 * Settings:
 * @param outputHistory	Should fire history be saved?  Defaults to false.
 *
 * @note The calculated output vectors do not need to be sized on input (though they can be).  Empty
 * vectors can be passed in, which simplifies the calling code.  The vectors will be resized on
 * return.
 *
 * @returns On return many arguments may be sorted, updated, or returned. [MORE!!!!!]
 *
 * @note The number of fuel types could be inferred from the input data eliminating the need for the
 *       number parameter.  The number parameter is currently used to confirm the other inputs are
 *       consistent.
 *
 * @par History:
 * Added as an programatic alternative entry point to the original interactive program.
 */
void Simulate(double& fi, const double ti, const double u, const double d, const double tpamb,
              const double ak, const double r0, const double dr, double& dt, const double wdf,
              const double dfm, const int ntimes, const int number,
              std::vector<std::string>& parts, std::vector<double>& wdry, std::vector<double>& ash,
              std::vector<double>& htval, std::vector<double>& fmois, std::vector<double>& dendry,
              std::vector<double>& sigma, std::vector<double>& cheat, std::vector<double>& condry,
              std::vector<double>& tpig, std::vector<double>& tchar, std::vector<double>& wo,
              std::vector<double>& xmat, std::vector<double>& tign, std::vector<double>& tout,
              std::vector<double>& diam, const bool outputHistory)
{
	int lenkl = Length_kl(number);//Triangular matrix size.

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
	std::vector<double> tdry(lenkl);		//Time of drying start of the larger of each fuel component pair, s. [maxkl]
	std::vector<double> wodot(lenkl);		//Dry loading loss rate for larger of pair. [maxkl]
	std::vector<double> ddot(lenkl);		//Diameter reduction rate, larger of pair, m / s. [maxkl]
	std::vector<double> qcum(lenkl);		//Cumulative heat input to larger of pair, J / sq m. [maxkl]
	std::vector<double> tcum(lenkl);		//Cumulative temp integral for qcum (drying). [maxkl]
	std::vector<double> acum(lenkl);		//Heat pulse area for historical rate averaging. [maxkl]
	std::vector<std::vector<double>> qdot(lenkl, std::vector<double>(mxstep));//History (post ignite) of heat transfer rate
	                                    	//to the larger of each component pair, W / sq m. [maxkl, mxstep]
	std::vector<int> key(number);			//Ordered index list. [maxno]
	std::vector<std::string> list(number);	//Intermediary for reordering parts name array. [maxno]
	                                 		//Probably not needed here.  See notes in ARRAYS().

	//Scalars in order of appearance:
	double dfi; 	//Duff fire intensity (aka I sub d) for DUFBRN().
	double tdf; 	//Burning duration (aka t sub d) for DUFBRN().
	int now;		//Index marking the current time step.
	double tis;		//Current time (ti + number of time steps * dt).
	int ncalls;		//Counter of calls to START().
	double fid;		//Fire intensity due to duff burning.

	//In the original code fmin was a local treated as a constant.  Passing it in might be good:
	const double fimin = 0.1;//Fire intensity (kW / sq m) at which fire goes out.

	//There are a large number of locals in this routine that are not explictly initialized.
	//Most are initialized in ARRAYS() and START().  Testing was done to confrim that explicit
	//initialization was not needed for the remaining locals.

	//Set SaveHistory:
	SaveHistory = outputHistory;

	//Validate the size of incoming input vectors:
	if (parts.size() != number || wdry.size() != number ||
	    !SameLengths(wdry, ash, htval, fmois, dendry, sigma, cheat) ||
	    !SameLengths(wdry, condry, tpig, tchar))
	{
		Stop("Fuel properties must be of equal length.");
	}
	NumFuelTypes = number;
	
	//Validate and correct the size of output vectors:
	ValidateOutputVector(wo, "wo");
	ValidateOutputVector(xmat, "xmat");
	ValidateOutputVector(tign, "tign");
	ValidateOutputVector(tout, "tout");
	ValidateOutputVector(diam, "diam");

	//Sort the fuel components and calculate the interaction matrix...
	ARRAYS(wdry, ash, dendry, fmois, sigma, htval, cheat, condry, tpig, tchar, diam, key, work, ak,
	       elam, alone, xmat, wo, parts, list, area);

	//Record the state before the start of the simulation.  This need to be done after ARRAYS()
	//because parts, wo, and diam may get reordered.  now and tis are not initialized yet so we
	//set the time explicitly.  
	SaveStateToFile(0, 0.0, number, parts, wo, diam, fi);
	//The first simulated time point starts at the end of the igniting fire residence time.  The
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

	if (tign[0] < 0.0)//Fuels failed to ignite. [Review: This is not completely consistent with the notes above?????]
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

	//Record the state after START() and the first call to FIRINT(), if needed:
	SaveStateToFile(now, tis, number, parts, wo, diam, fi);

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

			//Save the state at each timestep, if needed:
			SaveStateToFile(now, tis, number, parts, wo, diam, fi);

			if (fi <= fimin)
			{
				break;
			}
		}
	}//(fi > fimin)

	//Return the time when the fire dropped below fimin as the time the fire "went out".
	//We return the value in dt.  This may be changed in the future.
	dt = tis;
	//We could also return the number of timesteps completed in ntimes but that doesn't add much.
}

//SimulateR():
/** This is a wrapper for Simulate() that allows it to be called from R (using .C()):
 *
 * The parameters parallel those of Simulate() but are pointers meet the requirements of R's .C()
 * function.  There is no easy way to pass a string so the 'parts' parameter is omitted.
 *
 * Igniting fire and environmental data:
 * @param[in, out] fi	Current fire intensity (site avg), kW / sq m
 * @param[in] ti		Igniting fire residence time (s).
 * @param[in] u			Mean horizontal windspeed at top of fuelbed (m/s).
 * @param[in] d			Fuelbed depth (m).
 * @param[in] tpamb		Ambient temperature (K).
 * 
 * Internal and control variables:
 * @param[in] ak		Area influence factor (ak / K_a parameter).
 *              		We modify the original behavior such that a negative value indicates that
 *              		the value of ak / K_a should be calculated according to Albini & Reinhardt
 *              		1997.
 * @param[in] r0		Minimum value of mixing parameter.
 * @param[in] dr		Max - min value of mixing parameter.
 * @param[in] dt		Time step for integration of burning rates (s).
 *              		On completion contains the time the fire went out.
 *              		A value of -1 indicates the fuel did not ignite.  A value
 *              		of -2 indicates the fuel did not complete drying.  In such
 *              		cases most of remaining return variables will not be
 *              		meaningful.
 * 
 * @param[in] wdf		Duff loading (kg/m^2, aka W sub d).
 * @param[in] dfm		Ratio of moisture mass to dry organic mass /
 *               		duff fractional moisture (aka R sub M).
 * @param[in] ntimes	Number of time steps to run.  Move down?
 * @param[in] number	The number of fuel classes.
 *
 * Fuel component property arrays:
 * @param[in,out] wdry		Ovendry mass loading, kg/sq m. [maxno]
 * @param[in,out] ash		Mineral content, fraction dry mass. [maxno]
 * @param[in,out] htval		Low heat of combustion (AKA heat content), J / kg. [maxno]
 * @param[in,out] fmois		Moisture fraction of component. [maxno]
 * @param[in,out] dendry	Ovendry mass density, kg / cu m. [maxno]
 * @param[in,out] sigma		Surface to volume ratio, 1 / m. [maxno]
 * @param[in,out] cheat		Specific heat capacity, (J / K) / kg dry mass. [maxno]
 * @param[in,out] condry	Thermal conductivity, W / m K, ovendry. [maxno]
 * @param[in,out] tpig		Ignition temperature, K. [maxno]
 * @param[in,out] tchar		Char temperature, K. [maxno]
 *
 * Calculated outputs:
 * JMR_Note: No longer in argument order!!!!!
 * @param[out] wo		Current ovendry loading for the larger of each component pair, kg/sq m. [maxkl]
 * @param[out] xmat		Table of influence fractions between components. [maxkl]
 * @param[out] tign		Ignition time for the larger of each fuel component pair, s. [maxkl]
 * @param[out] tout		Burnout time of larger component of pairs, s. [maxkl]
 * @param[out] diam		Current diameter of the larger of each fuel component pair, m. [maxkl]
 *
 * @param[in] outputHistory	Should fire history be saved? (0 = no, 1 = yes)
 *
 * @par History:
 * Added for module.
 */
extern "C" void SimulateR(double* fi, const double* ti, const double* u, const double* d,
                          const double* tpamb, const double* ak, const double* r0, const double* dr,
                          double* dt, const double* wdf, const double* dfm, const int* ntimes,
                          const int* number, double* wdry, double* ash, double* htval,
                          double* fmois, double* dendry, double* sigma, double* cheat,
                          double* condry, double* tpig, double* tchar, double* wo, double* xmat,
                          double* tign, double* tout, double* diam, const int* outputHistory)
{
	//Local type conversion intermediates:
	bool historyLogical;
	
	//Character strings can't be passed in from R so we assemble some generic names to pass in:
	std::vector<std::string> parts(*number);//Fuel component names / labels. [maxno]

	for (int i = 0; i < *number; i++)
	{
		parts[i] = "Fuel " + std::to_string(i);
	}

	//Convert input arrays to vectors:
	std::vector<double> wdryVec(wdry, wdry + *number);
	std::vector<double> ashVec(ash, ash + *number);
	std::vector<double> htvalVec(htval, htval + *number);
	std::vector<double> fmoisVec(fmois, fmois + *number);
	std::vector<double> dendryVec(dendry, dendry + *number);
	std::vector<double> sigmaVec(sigma, sigma + *number);
	std::vector<double> cheatVec(cheat, cheat + *number);
	std::vector<double> condryVec(condry, condry + *number);
	std::vector<double> tpigVec(tpig, tpig + *number);
	std::vector<double> tcharVec(tchar, tchar + *number);

	//Triangular matrix size, equivalent to maxkl but determined by 'number' passed in:
	int lenkl = Length_kl(*number);
	
	std::vector<double> xmatVec(xmat, xmat + lenkl);
	std::vector<double> tignVec(tign, tign + lenkl);
	std::vector<double> toutVec(tout, tout + lenkl);
	std::vector<double> woVec(wo, wo + lenkl);
	std::vector<double> diamVec(diam, diam + lenkl);

	//R logical variables come in as ints:
	if (*outputHistory == 0)
	{
		historyLogical = false;
	}
	else
	{
		historyLogical = true;//The value should be 1.  We don't check for the NA value or others.
	}

	Simulate(*fi, *ti, *u, *d, *tpamb, *ak, *r0, *dr, *dt, *wdf, *dfm, *ntimes, *number,
	         parts,
	         wdryVec, ashVec, htvalVec, fmoisVec, dendryVec,
	         sigmaVec, cheatVec, condryVec, tpigVec, tcharVec,
	         woVec, xmatVec, tignVec, toutVec, diamVec,
	         historyLogical);

	//Convert outputs back arrays:
	for (int i = 0; i < *number; i++)
	{
		wdry[i] = wdryVec[i];
		ash[i] = ashVec[i];
		htval[i] = htvalVec[i];
		fmois[i] = fmoisVec[i];
		dendry[i] = dendryVec[i];
		sigma[i] = sigmaVec[i];
		cheat[i] = cheatVec[i];
		condry[i] = condryVec[i];
		tpig[i] = tpigVec[i];
		tchar[i] = tcharVec[i];
		xmat[i] = xmatVec[i];
		tign[i] = tignVec[i];
		tout[i] = toutVec[i];
		wo[i] = woVec[i];
		diam[i] = diamVec[i];
	}
	
	for (int j = 0; j < lenkl; j++)
	{
		xmat[j] = xmatVec[j];
		tign[j] = tignVec[j];
		tout[j] = toutVec[j];
		wo[j] = woVec[j];
		diam[j] = diamVec[j];
	}
}

//DUFBRN():
/** Calculate duff fire properties.
 *
 * @par Original Burnup Description:
 * !c Duff burning rate (ergo, intensity) and duration
 *
 * @param[in] wdf	Duff loading (kg/m^2, aka W sub d)
 * @param[in] dfm	Ratio of moisture mass to dry organic mass / duff fractional moisture (aka R sub M)
 * @param[out] dfi	Duff fire intensity (aka I sub d)		Units?????
 * @param[out] tdf	Burning duration (aka t sub d)			Units?????
 *
 * @returns On return dfi and tdf are returned in parameters.
 *
 * @par History:
 * Modernized original Burnup subroutine.
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

//DufBrnR() probably does not need to be ported.

//ARRAYS():
/** Prepare arrays.  ...
 *
 * @par Original Burnup Description:
 * !c Orders the fuel description arrays according to the paradigm described in
 * !c subroutine SORTER and computes the interaction matrix xmat from the array
 * !c elam and the list alone returned from subroutine OVLAPS.
 *
 * The large number of parameters are hard to handle...

 * @param[in,out] wdry		Ovendry mass loading, kg / sq m. [maxno]
 * @param[in,out] ash		Mineral content, fraction dry mass. [maxno]
 * @param[in,out] dendry	Ovendry mass density, kg / cu m. [maxno]
 * @param[in,out] fmois		Moisture content, fraction dry mass. [maxno]
 * @param[in,out] sigma		Surface to volume ratio, 1 / m. [maxno]
 * @param[in,out] htval		Low heat of combustion, J / kg. [maxno]
 * @param[in,out] cheat		Specific heat capacity, (J / K) / kg dry mass. [maxno]
 * @param[in,out] condry	Thermal conductivity, W / m K, ovendry. [maxno]
 * @param[in,out] tpig		Ignition temperature, K. [maxno]
 * @param[in,out] tchar		Char temperature, K. [maxno]
 * @param[out] diam			Initial diameter, m [by interaction pairs]. [maxkl]
 * @param[out] key 			Ordered index list. [maxno]
 * @param[out] work			Workspace array. [maxno]
 * @param[in]  ak			Area influence factor [ak parameter].
 * @param[out] elam			Interaction matrix from OVLAPS. [maxno, maxno]
 * @param[out] alone		Noninteraction fraction list from OVLAPS. [maxno]
 * @param[out] xmat			Consolidated interaction matrix. [maxkl]
 * @param[out] wo			Initial dry loading by interaction pairs. [maxkl]
 * @param[in,out] parts		Fuel component names / labels. [maxno]
 * @param[out] list			Intermediary for reordering parts name array. [maxno]
 *                 			This is passed in but is not initialized prior
 *                 			to that.  It doesn't appear that it is used
 *                 			after it is returned.  It appears to only be
 *                 			used internal to this routine.
 *                 			The same seems to be true for elam and alone?
 * @param[in,out] area		Fraction of site area expected to be covered at
 *                    		least once by initial planform area of ea size. [maxno]
 * @param[in] number		The actual number of fuel classes.  If omitted
 *                  		this will be determined from the other inputs.
 *
 * @returns On return many arguments may be sorted, updated, or returned. [MORE!!!!!]
 *
 * @par History:
 * Modernized original Burnup subroutine.
 * Several arguments have been removed that were present in the original routine.  The number
 * argument has been moved and is now optional and is only used in the interactive context.
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
	//The original Fortran code used counters k and j, which are 1 bases indexes.  We use k and j
	//for 1 based triangular matrix indexes and k0 and j0 for 0 based array indexes.
	int k0;//Counter

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

	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		k0 = key[j0];
		list[j0] = parts[k0];
	}
	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		parts[j0] = list[j0];
	}

	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		k0 = key[j0];
		work[j0] = wdry[k0];
	}
	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		wdry[j0] = work[j0];
	}

	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		k0 = key[j0];
		work[j0] = ash[k0];
	}
	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		ash[j0] = work[j0];
	}

	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		k0 = key[j0];
		work[j0] = htval[k0];
	}
	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		htval[j0] = work[j0];
	}

	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		k0 = key[j0];
		work[j0] = cheat[k0];
	}
	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		cheat[j0] = work[j0];
	}

	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		k0 = key[j0];
		work[j0] = condry[k0];
	}
	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		condry[j0] = work[j0];
	}

	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		k0 = key[j0];
		work[j0] = tpig[k0];
	}
	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		tpig[j0] =  work[j0];
	}

	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		k0 = key[j0];
		work[j0] = tchar[k0];
	}
	for (int j0 = 0; j0 < numFuelTypes; j0++)
	{
		tchar[j0] = work[j0];
	}

	OVLAPS(wdry, sigma, dendry, ak, fmois, xmat, elam, alone, area, number);

	for (int k = 1; k <= numFuelTypes; k++)
	{
		k0 = k - 1;
		double diak = 4.0 / sigma[k0];
		double wtk = wdry[k0];

		//Populate the alone/no companion indexes of the arrays:
		int kl = Loc(k, 0);
		diam[kl] = diak;
		xmat[kl] = alone[k0];
		wo[kl] = wtk * xmat[kl];

		//Populate the interacting indexes of the arrays:
		for (int j = 1; j <= k; j++)
		{
			int kj = Loc(k, j);
			diam[kj] = diak;
			xmat[kj] = elam[k0][j - 1];//Convert j to 0 based index, j0.
			wo[kj] = wtk * xmat[kj];
		}
	}
}

//SORTER():
/** Sort the fuels.
 *
 * @par Original Burnup Description:
 * !c Sorts fuel element list in order of increasing size (decreasing sigma)
 * !c For elements with same size order determined on increasing moisture
 * !c content (fmois). If items have same size and moisture content, order
 * !c on the basis of increasing mass density (dryden). "number" elements are
 * !c included in the list, which has a maximum length of "maxno". The integer
 * !c list: key(j), j = 1, number holds the indices in order, so other
 * !c fuel parameters can be ordered and associated as necessary.
 *
 * @param[in,out] sigma(:)		Surface to volume ratio, 1 / m. [maxno]
 * @param[in,out] fmois(:)		Moisture content, fraction dry mass. [maxno]
 * @param[in,out] dryden(:)		Ovendry mass density, kg / cu m. [maxno]
 * @param[in,out] key(:)		Ordered index list. [maxno]					//out rather than inout?????
 * @param[in] number			The actual number of fuel classes.  If omitted
 *                  			this will be determined from the other inputs.
 * Since SORTER() is only downstream of ARRAYS() making number optional doesn't gain us much.
 *
 * @returns On return sigma, fmois, dryden are reordered and key is set.
 *
 * @par History:
 * Modernized original Burnup subroutine.
 * The maxno argument has been removed.  The number argument has been moved and is now optional.
 */
void SORTER(std::vector<double>& sigma, std::vector<double>& fmois, std::vector<double>& dryden,
            std::vector<int>& key, const int number)
{
	int maxNumFuelTypes;		//The maximum number of fuel classes allowed. The input arrays
								//may exceed the actual number. (Not present in original code.)
	int numFuelTypes;			//The actual number of fuel types, explicit or implied.
	int i0;						//Counter.
	double s, fm, de, keep, usi;//Hold the values of the current search index.  s and usi = inverse of SAV.
	bool diam, mois, dens, tied;
	bool newIndexFound;//Note: Not present in original code.

	maxNumFuelTypes = sigma.size();//Determine maximum number of fuels from input arrays.
	
	//Determine the actual number of fuel types:
	if (number > 0)
	{
		numFuelTypes = number;
	}
	else
	{
		numFuelTypes = maxNumFuelTypes;
	}

	newIndexFound = false;

	for (int j0 = 0; j0 < maxNumFuelTypes; j0++)
	{
		key[j0] = j0;
	}

	//!c Replacement sort: order on increasing size, moisture, density
	for (int j0 = 1; j0 < numFuelTypes; j0++)
	{
		//Store the values for this fuel index:
		s = 1.0 / sigma[j0];
		fm = fmois[j0];
		de = dryden[j0];
		keep = key[j0];

		//Compare this index (j0) with every index before it:
		for (i0 = (j0 - 1); i0 >= 0; i0--)//i is only used as a 0 based array index.
		{
			usi = 1.0 / sigma[i0];
			diam = (usi < s);

			if (diam)
			{
				newIndexFound = true;
				break;
			}

			tied = (usi == s);

			if (tied)
			{
				mois = (fmois[i0] < fm);
				if (mois)
				{
					newIndexFound = true;
					break;
				}

				tied = (fmois[i0] == fm);
				if (tied)
				{
					dens = (dryden[i0] <= de);
					if (dens)
					{
						newIndexFound = true;
						break;
					}
				}
			}

			//i is greater than j.
			//Move entry i (the entry we are comparing to) down one:
			sigma[i0 + 1] = sigma[i0];
			fmois[i0 + 1] = fmois[i0];
			dryden[i0 + 1] = dryden[i0];
			key[i0 + 1] = key[i0];
		}

		if (newIndexFound)
		{
			//If a new location has been identified record entry j at i + 1 (below).
			newIndexFound = false;//Reinitialize for the next comparison.
		}
		else
		{
			i0 = -1;//Reseting i to -1 will move this entry (j) to postion 0.
		}

		//Record the values for fuel index j0 at the identified index:
		sigma[i0 + 1] = 1.0 / s;
		fmois[i0 + 1] = fm;
		dryden[i0 + 1] = de;
		key[i0 + 1] = keep;
	}
}

//OVLAPS():
/** Compute the fuel interaction matrix.
 *
 * @par Original Burnup Description:
 * !c Computes the interaction matrix elam(j, k) which apportions the
 * !c influence of smaller and equal size pieces on eaqh size class for the
 * !c purpose of establishing the rates at which the elemnts burn out.
 * !c Input quantities are dryld, the ovendry mass per unit area of each
 * !c element available for burning after the passage of the igniting surface
 * !c fire; sigma, the particle's surface/ volume ratio, and dryden, the
 * !c ovendry mass density of the particle; ak a dimensionless parameter that
 * !c scales the planform area of a particle to its area of influence. There
 * !c are "number" separate particle classes, of a maximum number = maxno.
 * !c It is assumed that the 1ist are ordered on size class (nonincreasing
 * !c surface/ volume ratio). List "alone" gives the fraction of each loading
 * !c that is not influenced by any other category.
 *
 * @param[in] dryld		Ovendry mass per unit area of each element (kg/sq m) (= wdry, ...). [maxno]
 * @param[in] sigma		Surface to volume ratio, 1 / m. [maxno]
 * @param[in] dryden	Ovendry mass density, kg / cu m (elsewhere dendry). [maxno]
 * @param[in] ak		Area influence factor (ak / K_a parameter).
 * @param[in] fmois		Moisture fraction of component. [maxno]
 * @param[out] beta		Consolidated interaction matrix. (elsewhere = xmat). [maxkl]
 * @param[out] elam		Interaction matrix. [maxno, maxno]
 * @param[out] alone	Non-interacting fraction for each fuel class. [maxno]
 * @param[out] area		Fraction of site area expected to be covered at
 *                 		least once by initial planform area of ea size. [maxno]
 * @param[in] number	The actual number of fuel classes.  If omitted
 *                  	this will be determined from the other inputs.
 * Since OVLAPS() is only downstream of ARRAYS() making number optional doesn't gain us much.
 *
 * @returns On return beta, elam, alone, and area are returned in parameters.
 *
 * @par History:
 * Modernized original Burnup subroutine.
 * We modify the original behavior such that a negative value for ak indicates that the value of
 * ak / K_a should be calculated.  This requires fmois to be passed in, which was not one of the
 * original arguments.
 * Several arguments have been removed that were present in the original routine.  The number
 * argument has been moved and is now optional. 
 */
void OVLAPS(const std::vector<double> dryld, const std::vector<double> sigma,
            const std::vector<double> dryden, const double ak, const std::vector<double> fmois,
            std::vector<double>& beta, std::vector<std::vector<double>>& elam,
            std::vector<double>& alone, std::vector<double>& area, const int number)
{
	double pi;//Convert to a constant?
	int numFuelTypes;//The actual number of fuel types, explicit or implied.
	int kl;//Counter
	double siga;//K_a * diameter
	double a;
	double bb;
	double frac;//Intermediate for calculating alone().
	double K_a;//Value of K_a (ak) parameter, fixed or calculated depending on the mode.

	pi = std::abs(std::acos(-1.0));//Calculate pi.

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
	std::fill(alone.begin(), alone.end(), 0);
	std::fill(beta.begin(), beta.end(), 0);
	for (auto& row : elam)
	{
		std::fill(row.begin(), row.end(), 0);
	}
	//We should allow the vectors of any size, including empty vectors, to be passed in!!!!!

	for (int k = 1; k <= numFuelTypes; k++)
	{
		int k0 = k - 1;

		for (int l = 1; l <= k; l++)//Use 0 indexing.
		{
			int l0 = l - 1;

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
				//K_a = 3.25 * exp(-20 * pow(fmois[l], 2));
				K_a = 3.25 * exp(-20 * pow(fmois[l0], 2));
			}

			//SAV / pi = diameter (units are carried by ak):
			siga = K_a * sigma[k0] / pi;

			kl = Loc(k, l);
			a = siga * dryld[l0] / dryden[l0];//siga * ? units in meters
			if (k == l)
			{
				bb = 1.0 - exp(-a);			// JMR: FOFEM suggests this can hit 0?
				area[k0] = bb;
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
		elam[0][0] = beta[1];
		alone[0] = 1.0 - elam[0][0];
		return;
	}

	for (int k = 1; k <= numFuelTypes; k++)
	{
		//These inner loops could be combined to simplify the logic and make it more readable!!!!!
		int k0 = k - 1;
		
		frac = 0.0;
		for (int l = 1; l <= k; l++)
		{
			kl = Loc(k, l);
			frac = frac + beta[kl];
		}

		if (frac > 1.0)
		{
			for (int l = 1; l <= k; l++)
			{
				kl = Loc(k, l);
				elam[k0][l - 1] = beta[kl] / frac;//l0
			}
			alone[k0] = 0.0;
		}
		else
		{
			for (int l = 1; l <= k; l++)
			{
				kl = Loc(k, l);
				elam[k0][l - 1] = beta[kl];//l0
			}
			alone[k0] = 1.0 - frac;
		}
	}
}

//START():
/** Start the fuel consumption simulation.
 *
 * @par Original Burnup Description:
 * !c This routine initializes variables prior to starting sequence of calls
 * !c to subroutine STEP.  On input here, fi is area intensity of spreading
 * !c fire, dt is the residence time for the spreading fire.  Only physical
 * !c parameters specified are the fuel array descriptors. To call STEP,
 * !c one must initialize the following variables.
 *
 * Arguments: (by category, not argument order)

! JMR_NOTE: The original comments imply that alfa, diam, and wo should all be intent(in).
! However the code is not consistant with that.

 * @param[in] dt		Spreading fire residence time (s) (= ti, tis, or time elsewhere).
 * @param[in] now		Index marks end of time step.
 * @param[in,out] wo	Current ovendry loading for the larger of each component pair, kg / sq m.
 *                  	Updated on return. [maxkl]
 * @param[out] alfa		Dry thermal diffusivity of component, sq m / s. [maxno]
 * @param[in] dendry	Ovendry density of component, kg / cu m. [maxno]
 * @param[in] fmois		Moisture fraction of component. [maxno]
 * @param[in] cheat		Specific heat capacity of component, J / kg K. [maxno]
 * @param[in] condry	Ovendry thermal conductivity, W / sq m K. [maxno]
 * @param[in,out] diam	Current diameter of the larger of each
 *                    	fuel component pair, m.  Updated on return. [maxkl]
 * @param[in] tpig		Ignition temperature (K), by component. [maxno]
 * @param[in] tchar		tchar = end - pyrolysis temperature (K), by component. [maxno]
 * @param[in] xmat		Table of influence fractions between components. [maxkl]
 * @param[in] tpamb		Ambient temperature (K).
 * @param[in] fi		Current fire intensity (site avg), kW / sq m.
 *
 * Parameters updated (input and output):
 * @param[in,out] ncalls	Counter of calls to this routine:
 *                      	= 0 on first call or reset,
 *                      	cumulates after first call.
 *                      	JMR_NOTE: This is a strange argument as it is only
 *                      	initialized to zero here and is not used.  It is
 *                      	returned and passed on to STEP().  It is probably
 *                      	initialized here to prevent isses related to
 *                      	persistance should more than one simulation be run in
 *                      	an interactive session.
 * @param[out] flit		Fraction of each component currently alight. [maxno]
 * @param[out] fout		Fraction of each component currently gone out. [maxno]
 * @param[out] tdry		Time of drying start of the larger of each fuel component pair, s. [maxkl]	[Units = s?????]
 * @param[out] tign		Ignition time for the larger of each fuel component pair, s. [maxkl]
 * @param[out] tout		Burnout time of larger component of pairs, s. [maxkl]
 * @param[out] qcum		Cumulative heat input to larger of pair, J / sq m. [maxkl]
 * @param[out] tcum		Cumulative temp integral for qcum (drying). [maxkl]
 * @param[out] acum		Heat pulse area for historical rate averaging. [maxkl]
 * @param[out] qdot		History (post ignite) of heat transfer rate
 *                 		to the larger of each component pair, W / sq m. [maxkl, mxstep]
 * @param[out] ddot		Diameter reduction rate, larger of pair, m / s. [maxkl]
 * @param[out] wodot	Dry loading loss rate for larger of pair. [maxkl]
 * @param[in,out] work	Workspace array. [maxno]		[Or alternative description in STEP()?????]
 *
 * Constant parameters:
 * @param[in] u			Mean horizontal windspeed at top of fuelbed (m/s).
 * @param[in] d			Fuelbed depth (m).
 * @param[in] r0		Minimum value of mixing parameter.
 * @param[in] dr		Max - min value of mixing parameter.
 *
 * Interactive context:
 * @param[in] number	The actual number of fuel classes.  If omitted
 *                  	this will be determined from the other inputs.
 *                  	Here the presence of number is used to determine
 *                  	if we are in an interactive session as well.
 *
 * @returns Nothing. ...
 
 * @par History:
 * Modernized original Burnup subroutine.
 * Several arguments have been removed that were present in the original routine.  The number
 * argument has been moved and is now optional and is only used in the interactive context.
 *
 * @note The constants ch2o and tpdry were included as arguments in the original code.  They have
 * chnaged to globals.
 * @note The original comments include hvap as a constant, but is not actually used:
 * hvap = heat of vaporization of water J / kg
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

	for (int k = 1; k <= numFuelTypes; k++)
	{
		int k0 = k - 1;

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

		for (int l = 0; l <= k; l++)//l in kl space, 0 based
		{
			kl = Loc(k, l);
			tout[kl] = rindef;
			tign[kl] = rindef;
			tdry[kl] = rindef;
			tcum[kl] = 0.0;
			qcum[kl] = 0.0;
		}
	}

	//The original code does not initialize acum and qdot.  Failure to do so leads to small
	//variations between runs and changes to results.
	std::fill(acum.begin(), acum.end(), 0);
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
	if (tf <= (tpdry + 10.0))
	{
		if (number > 0)//In an interactive session, preserve the original behavior:
		{
			Stop("Igniting fire cannot dry fuel");
		}
		else//Otherwise signal the condition and return
		{
			std::fill(tign.begin(), tign.end(), -2.0);//Value signals fuel did not dry.
			Msg.Log("Igniting fire cannot dry fuel.");
			return;
		}
	}
	thd = (tpdry - ts) / (tf - ts);
	tx = 0.5 * (ts + tpdry);

	for (int k = 1; k <= numFuelTypes; k++)//k0 = k as base 0 index
	{
		int k0 = k - 1;

		factor = dendry[k0] * fmois[k0];
		conwet = condry[k0] + 4.27e-04 * factor;
		for (int l = 0; l <= k; l++)//l in kl space, 0 based
		{
			kl = Loc(k, l);
			dia = diam[kl];
			HEATX(u, d, dia, tf, tx, hf, hb, conwet, en);
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

	for (int k = 1; k <= numFuelTypes; k++)
	{
		int k0 = k - 1;

		c = condry[k0];
		tigk = tpig[k0];
		for (int l = 0; l <= k; l++)
		{
			kl = Loc(k, l);
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

	for (int k = 1; k <= numFuelTypes; k++)
	{
		int k0 = k - 1;//Only used once!!!!!

		if (flit[k0] > 0.0)
		{
			nlit = nlit + 1;
		}

		for (int l = 0; l <= k; l++)
		{
			kl = Loc(k, l);
			trt = std::min(trt, tign[kl]);
		}
	}

	if (nlit == 0)
	{
		if (number > 0)//In an interactive session, preserve the original behavior:
		{
			Stop("START ignites no fuel");
		}
		else//Otherwise signal the condition and return:
		{
			std::fill(tign.begin(), tign.end(), -1.0);//Value signals fuel did not ignite.
			Msg.Log("Igniting fire cannot ignite fuel.");
			return;
		}
	}

	//!c Deduct trt from all time estimates, resetting time origin

	for (int k = 1; k <= numFuelTypes; k++)
	{
		for (int l = 0; l <= k; l++)
		{
			kl = Loc(k, l);
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

	for (int k = 1; k <= numFuelTypes; k++)
	{
		int k0 = k - 1;

		if (flit[k0] == 0.0)
		{
			for (int l = 0; l <= k; l++)
			{
				kl = Loc(k, l);
				ddot[kl] = 0.0;
				tout[kl] = rindef;
				wodot[kl] = 0.0;
			}
		}
		else
		{
			ts = tchar[k0];
			c = condry[k0];
			for (int l = 0; l <= k; l++)
			{
				kl = Loc(k, l);
				dia = diam[kl];
				HEATX(u, d, dia, tf, ts, hf, hb, c, e);
				qdot[kl][now - 1] = hb * std::max((tf - ts), 0.0);//Convert to 0 based index.
				aint = pow((c / hb), 2);
				ddt = dt - tign[kl];
				acum[kl] = aint * ddt;
				ddot[kl] = qdot[kl][now - 1] * work[k0];//Convert to 0 based index.
				tout[kl] = dia / ddot[kl];
				dnext = std::max(0.0, (dia - ddt * ddot[kl]));
				wnext = wo[kl] * pow((dnext / dia), 2);
				wodot[kl] = (wo[kl] - wnext) / ddt;

				diam[kl] = dnext;
				wo[kl] = wnext;
				df = 0.0;
				if (dnext <= 0.0)
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

//FIRINT():
/** Compute fire intensity.
 *
 * @par Original Burnup Description:
 * !c Computes fi = site avg fire intensity given the burning rates of all
 * !c interacting pairs of fuel components [ wodot ], the mineral ash content
 * !c of each component [ ash ], the heat of combustion value [ htval ] for
 * !c each, and the number of fuel components [ number ], where max - maxno.
 * !c fi is in kW / sq m, while htval is in J / kg.
 *
 * !c fint(k) is the correction to fi to adjust
 * !c the intensity level to be the local value where size k is burning.
 *
 * @param[in] wodot		Burning rates of interacting pairs of fuel components. [maxkl]
 * @param[in] ash		Mineral content, fraction dry mass. [maxno]
 * @param[in] htval		Low heat of combustion, J / kg. [maxno]
 * @param[in] number	The actual number of fuel classes.
 * @param[in] area		Fraction of site area expected to be covered at
 *                		least once by initial planform area of ea size. [maxno]
 * @param[out] fint		A vector to hold the corrected local fire intensity for each fuel type. [maxno]
 * @param[out] fi		Site avg fire intensity (kW / sq m).
 *
 * @returns fint[] and fi are returned it the parameters.
 *
 * @par History:
 * Modernized original Burnup subroutine.
 * Several arguments have been removed that were present in the original routine.
 */
void FIRINT(const std::vector<double> wodot, const std::vector<double> ash,
            const std::vector<double> htval, const int number, const std::vector<double> area,
            std::vector<double>& fint, double& fi)
{
	double sum;//Running total for fi.
	double wdotk;
	double term;
	double ark;//Area for element k.

	//Constants:
	const double small = 1.0e-06;

	//We could check that the inputs are the same size?????

	//Don't assume that fint is the right size!!!!!

	sum = 0.0;
	for (int k = 1; k <= number; k++)
	{
		int k0 = k - 1;

		wdotk = 0.0;
		for (int l = 0; l <= k; l++)
		{
			int kl = Loc(k, l);
			wdotk = wdotk + wodot[kl];
		}
		term = (1.0 - ash[k0]) * htval[k0] * wdotk * 1.0e-03;
		ark = area[k0];
		if (ark > small)
		{
			fint[k0] = term / ark - term;
		}
		else
		{
			fint[k0] = 0.0;
		}
		sum = sum + term;
	}

	fi = sum;
}

//TIGNIT():
/** This routine computes the halfspace surface ignition time under steady radiant heating with
 * surface film cooling.
 *
 * @param[in] tpam	Ambient temperature, K.
 * @param[in] tpdr	Fuel temperature at start of drying, K.
 *                	Currently this is always tpdry, so this argument could be cut.
 * @param[in] tpig	Fuel surface temperature at ignition, K.
 * @param[in] tpfi	Fire enviroriment temperature, K.
 * @param[in] cond	Fuel ovendry thermal conductivity, W / m K.
 * @param[in] chtd	Fuel ovendry specific heat capacity, J / kg K.
 * @param[in] fmof	Fuel moisture content, fraction dry weight.
 * @param[in] dend	Fuel ovendry density, kg / cu m.
 * @param[in] hbar	Effective film heat transfer coefficient [< HEATX] W / sq m K.
 *
 * @returns tmig Predicted time to piloted ignition, s.
 * 
 * @par History:
 * Modernized original Burnup subroutine.
 *
 * @note The orignal Fortran routine returns tmig as an argument.  !!!!!
 */
double TIGNIT(const double tpam, const double tpdr, const double tpig, const double tpfi,
              const double cond, const double chtd, const double fmof, const double dend,
              const double hbar)
{
	double b03;
	double xlo, xhi, xav;	//Binary search bounds and middle (average).
	double fav;				//Value of ff() for current search value.
	double beta, conw;
	double dtb;				//The temperature increase required to reach the drying temperature.
	double dti;				//The temperature increase required to reach the ignition temperature.
	double ratio, rhoc;
	double tmig;			//Return value: Predicted time to piloted ignition, s.

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
	// Testing was done to confrim that explicit initialization of other locals was not needed. [in Fortran]

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

		if (std::abs(fav) <= small)
		{
			break;
		}
		else if (fav < 0.0)
		{
			xlo = xav;
		}
		else if (fav > 0.0)
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
	tmig = (pow((beta / hbar), 2)) * conw * rhoc;

	return tmig;
}

//DRYTIM():
/** Compute the (dimensionless) drying time.
 *
 * @par Original Burnup Description:
 * !c Given a Nusselt number (enu, actually Biot number = h D / k)
 * !c and the dimensionless temperature rise required for the start
 * !c of surface drying (theta), returns the dimerisionless time (tau)
 * !c needed to achieve it. The time is given multiplied by thermal
 * !c diffusivity and divided by radius squared. Solution by binary search.
 *
 * @param[in] enu	Biot number.
 * @param[in] theta	Temperature rise required for the start of moisture loss.
 *
 * @returns tau The time required for the start of moisture loss.
 *
 * @par History:
 * Modernized original Burnup subroutine.
 *
 * @note The original Fortran routine returned tau via an argument.
 */
double DRYTIM(const double enu, const double theta)
{
	double xl, xh, xm;//The binary search low and high bounds, and the search center.
	double x;
	double approx;
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
	tau = pow((0.5 * x / enu), 2);

	return tau;
}


//HEATX():
/** This routine calculates how heat is transfered from the fire environment to a given fuel type.
 *
 * @par Original Burnup Description:
 * !c Given horizontal windspeed u at height d [top of fuelbed], cylindrical
 * !c fuel particle diameter dia, fire environment temperature tf, and mean
 * !c surface temperature, ts, subroutine returns film heat transfer coefficient
 * !c hfm and an "effective" film heat transfer coefficient including radiation
 * !c heat transfer, hbar.  Using the wood's thermal conductivity, cond, the
 * !c modified Nusselt number [ en ] used to estimate onset of surface drying
 * !c is returned as well.
 *
 * @param[in] u		Mean horizontal windspeed at top of fuelbed (m/s).
 * @param[in] d		Fuelbed depth (m).
 * @param[in] dia	Fuel diameter.									Units?????
 * @param[in] tf	Fire environment temperature.					Units?????
 * @param[in] ts	Mean surface temperature.						Units?????
 * @param[out] hfm	Film heat transfer coefficient.
 * @param[out] hbar	"Effective" film heat transfer coefficient.
 * @param[in] cond	Wood thermal conductivity.						Units?????
 * @param[out] en	Modified Nusselt number.
 *
 * @returns hfm, hbar, and en are returned via parameters.
 *
 * @par History:
 * Modernized original Burnup subroutine.
 *
 * @note The arguments are in the original order but it might be good to reorder so the return values are together!!!!!
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
		v = std::sqrt(u * u + 0.53 * g * d);
		re = v * dia / vis;
		enuair = 0.344 * pow(re, 0.56);
		conair = a + b * tf;
		fac = std::sqrt(std::abs(tf - ts) / dia);
		hfmin = fmfac * std::sqrt(fac);
		hfm = std::max((enuair * conair / dia), hfmin);
	}

	hrad = hradf * rad * (tf + ts) * (tf * tf + ts * ts);
	hbar = hfm + hrad;
	en = hbar * dia / cond;
}

//TEMPF():
/** Compute the fire environment temperature.
 *
 * !c Returns a fire environment temperature, TEMPF, given the fire intensity
 * !c q in kW / square meter, the ambient temperature tamb in Kelvins, and the
 * !c dimensionless mixing parameter r.
 *
 * @param[in] q Fire intensity (kW/m^2).
 * @param[in] r Dimensionless mixing parameter.
 * @param[in] tamb Ambient temperature (K).
 *
 * @returns Fire environment temperature Units????? (K)?????
 *
 * @par History:
 * Modernized original Burnup subroutine.
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

	while (true)
	{
		den = 1.0 + term * (rlast + 1.0) * (rlast * rlast + 1.0);
		rnext = 0.5 * (rlast + 1.0 + r / den);
		if (std::abs(rnext - rlast) < err)
		{
			tempf = rnext * tamb;
			break;
		}
		rlast = rnext;
	}

	return tempf;
}

//STEP():
/** This routine calculates one timestep of the fuel consumption process.
 *
 * @par Original Burnup Description:
 * !c Updates status of all fuel component pairs and returns a snapshot
 *
 * Arguments: (by category, not argument order)
 * @param[in] dt			Time step, sec.		[Or: Spreading fire residence time????? Confirm!!!!!]
 * @param[in] now			Index marks end of time step.
 * @param[in,out] wo		Current ovendry loading for the larger of each component pair, kg / sq m. [maxkl]
 * @param[in] alfa			Dry thermal diffusivity of component, sq m / s. [maxno]
 * @param[in] dendry		Ovendry density of component, kg / cu m. [maxno]
 * @param[in] fmois			Moisture fraction of component. [maxno]
 * @param[in] cheat			Specific heat capacity of component, J / kg K. [maxno]
 * @param[in] condry		Ovendry thermal conductivity, W / sq m K. [maxno]
 * @param[in,out] diam		Current diameter of the larger of each
 *                    		fuel component pair, m.  Updated on return. [maxkl]
 * @param[in] tpig			Ignition temperature (K), by component. [maxno]
 * @param[in] tchar			tchar = end - pyrolysis temperature (K), by component. [maxno]
 * @param[in] xmat			Table of influence fractions between components. [maxkl]
 * @param[in] tpamb			Ambient temperature (K).
 * @param[in] fi			Current fire intensity (site avg), kW / sq m.
 * @param[in] work			Factor of heat transfer rate hbar * (Tfire - Tebar)
 *                			that yields ddot (k). [maxno]
 *
 * Parameters updated (input and output):
 * @param[in,out] ncalls	Counter of calls to this routine:
 *                      	= 0 on first call or reset,
 *                      	cumulates after first call.
 * @param[in,out] flit		Fraction of each component currently alight. [maxno]
 * @param[in,out] fout		Fraction of each component currently gone out. [maxno]
 * @param[in,out] tdry		Time of drying start of the larger of each fuel component pair, s. [maxkl]
 * @param[in,out] tign		Ignition time for the larger of each fuel component pair, s. [maxkl]
 * @param[in,out] tout		Burnout time of larger component of pairs, s. [maxkl]
 * @param[in,out] qcum		Cumulative heat input to larger of pair, J / sq m. [maxkl]
 * @param[in,out] tcum		Cumulative temp integral for qcum (drying). [maxkl]
 * @param[in,out] acum		Heat pulse area for historical rate averaging. [maxkl]
 * @param[in,out] qdot		History (post ignite) of heat transfer rate
 *                    		to the larger of each component pair, W / sq m. [maxkl, mxstep]
 * @param[in,out] ddot		Diameter reduction rate, larger of pair, m / s. [maxkl]
 * @param[in,out] wodot		Dry loading loss rate for larger of pair. [maxkl]
 *
 * Constant parameters:
 * @param[in] u				Mean horizontal windspeed at top of fuelbed (m/s).
 * @param[in] d				Fuelbed depth (m).
 * @param[in] r0			Minimum value of mixing parameter.
 * @param[in] dr			Max - min value of mixing parameter.
 *
 * @param[in] tin			Start of current time step.
 * @param[in] fint			Correction to fi to compute local intensity
 *                			that may be different due to k burning. [maxno]
 * @param[in] fid			Fire intensity due to duff burning ... this is
 *               			used to up the fire intensity for fuel pieces
 *               			that are burning without interacting with others.
 *
 * Interactive context:
 * @param[in] number		The actual number of fuel classes.  If omitted
 *                  		this will be determined from the other inputs.

* Differences from the interface of START:
- alfa is input rather than output
- work is nput rather than output
- tin, fint, fid are added.

 * @par History:
 * Modernized original Burnup subroutine.
 * Several arguments have been removed that were present in the original routine.  The number
 * argument has been moved and is now optional and is only used in the interactive context.
 *
 * @note This routine takes a large number of arguments and the order is a bit confusing
 * with input and output parameters mixed in the order.
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
	int tIndex;		//Time index. = index in the Fortran code.  Renamed to avoid potential conflict with a POSIX function name.
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

	int kl;			//Triangular matrix index, 0 based.

	//There are a large number of locals in this routine that are not explictly initialized.
	//Testing was done to confrim that explicit initialization was not needed. [in Fortran]

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
	tifi = tnow - static_cast<double>(now - 1) * dt;
	next = now + 1;

	for (int k = 1; k <= numFuelTypes; k++)//Start major k loop.
	{
		int k0 = k - 1;

		c = condry[k0];
		for (int l = 0; l <= k; l++)//Start major l loop.
		{
			int l0 = l - 1;
			kl = Loc(k, l);
			tdun = tout[kl];

			//!c See if k of (k, l) pair burned out

			if (tnow >= tdun)
			{
				ddot[kl] = 0.0;
				wodot[kl] = 0.0;
				continue;//Back to start of major l loop.
			}
			if (tnext >= tdun)
			{
				tgo = tdun - tnow;
				ddot[kl] = diam[kl] / tgo;

				wodot[kl] = wo[kl] / tgo;
				wo[kl] = 0.0;
				diam[kl] = 0.0;
				continue;//Back to start of major l loop.
			}

			//!c k has not yet burned out ... see if k of (k, l) pair is ignited

			tlit = tign[kl];
			if (tnow >= tlit)
			{
				ts = tchar[k0];
				
				//In the original code the following conditionals were in series:
				if (l == 0)
				{
					r = r0 + 0.5 * dr;
					gi = fi + fid;
				}
				else if ((l != 0) && (l != k))//Or (l > 0) && (l < k).
				{
					r = r0 + 0.5 * (1.0 + flit[l0]) * dr;
					gi = fi + fint[k0] + flit[l0] * fint[l0];
				}
				else if (l == k)//Or just else.
				{
					r = r0 + 0.5 * (1.0 + flit[k0]) * dr;
					gi = fi + flit [k0] * fint[k0];
				}

				tf = TEMPF(gi, r, tpamb);
				dia = diam[kl];
				HEATX(u, d, dia, tf, ts, hf, hb, c, e);
				qqq = hb * std::max((tf - ts), 0.0);
				tst = std::max(tlit, tifi);
				
				nspan = std::max(l, static_cast<int>(std::round((tnext - tst) / dt)));//Ugly!!!!!
				if (nspan <= mxstep)
				{
					qdot[kl][nspan - 1] = qqq;//nspan is 1 based and must be converted.
				}
				else//if (nspan > mxstep)
				{
					for (int mu = 1; mu < mxstep; mu++)//mu has been adjusted to a 0 based index range.
					{
						qdot[kl][mu - 1] = qdot[kl][mu];
					}
					qdot[kl][mxstep - 1] = qqq;
				}
				aint = pow((c / hb), 2);
				acum[kl] = acum[kl] + aint * dt;

				//Time over which to perform averaging:
				tav1 = tnext - tlit;//Time since ignition.
				tav2 = acum[kl] / alfa[k0];//Measure of square of distance heat has penetrated fuel.
				tav3 = pow((dia / 4.0), 2) / alfa [k0];//Measure of time heat takes to reach center of fuel.
				tavg = std::min({tav1, tav2, tav3});

				tIndex = 1 + std::min(nspan, mxstep);
				qdsum = 0.0;
				tspan = 0.0;
				deltim = dt;

				//Calculate qdsum (sum of heat transfer (W/m^2 * s = J/m^2)):
				while (true)
				{
					tIndex = tIndex - 1;
					if (tIndex == 1)
					{
						deltim = tnext - tspan - tlit;
					}

					if ((tspan + deltim) >= tavg)
					{
						deltim = tavg - tspan;
					}

					qdsum = qdsum + qdot[kl][tIndex - 1] * deltim;//Convert to 0 index.
					tspan = tspan + deltim;

					if ((tspan < tavg) && (tIndex > 1))
					{
						continue;
					}
					else
					{
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
				continue;//Back to start of major l loop.
			}//if (tnow >= tlit)

			//!c See if k of (k, l) has reached outer surface drying stage yet

			dryt = tdry[kl];
			if ((tnow >= dryt) && (tnow < tlit))
			{
				//In the original code the following conditionals were in series:
				if (l == 0)
				{
					r = r0;
					gi = fi + fid;
				}
				if (l == k)
				{
					r = r0;
					gi = fi;
				}
				if ((l != 0) && (l != k))//Or just else.
				{
					r = r0 + 0.5 * flit[l0] * dr;
					gi = fi + flit[l0] * fint[l0];
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
					dtlite = TIGNIT(tpamb, tpdry, tpig[k0], tfe, condry[k0], cheat[k0], fmois[k0],
					                dendry[k0], heff);
				}
				tign[kl] = 0.5 * (dryt + dtlite);

				//!c If k will ignite before time step over, must interpolate

				if (tnext > tign[kl])
				{
					ts = tchar[k0];
					HEATX(u, d, dia, tf, ts, hf, hb, c, e);

					qdot[kl][0] = hb * std::max((tf - ts), 0.0);
					qd = qdot[kl][0];
					ddot[kl] = qd * work[k0];
					delt = tnext - tign[kl];
					dnext = std::max(0.0, dia - delt * ddot[kl]);
					wnext = wo[kl] * pow((dnext / dia), 2);
					if (dnext == 0.0)
					{
						tout[kl] = tnow + dia / ddot[kl];
					}
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
				continue;//Back to start of major l loop.
			}

			//!c If k of (k, l) still coming up to drying temperature, accumulate
			//!c heat input and driving temperature difference, predict drying start

			if (tnow < dryt)
			{
				factor = fmois[k0] * dendry[k0];
				conwet = condry[k0] + 4.27e-04 * factor;

				//In the original code the following conditionals were in series:
				if (l == 0)
				{
					r = r0;
					gi = fi + fid;
				}
				else if (l == k)
				{
					r = r0;
					gi = fi;
				}
				else if ((l != 0) && (l != k))//Or just else.
				{
					r = r0 + 0.5 * flit[l0] * dr;
					gi = fi + flit[l0] * fint[l0];
				}

				tf = TEMPF(gi, r, tpamb);
				if (tf <= (tpdry + 10.0))
				{
					continue;//Back to start of major l loop.
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
					continue;//Back to start of major l loop.
				}
				biot = he * dia / conwet;
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
						continue;//Back to start of major l loop.
					}

					dtlite = TIGNIT(tpamb, tpdry, tpig[k0], tf, condry[k0], cheat[k0], fmois[k0],
					                dendry[k0], hb);
					tign[kl] = 0.5 * (tdry[kl] + dtlite);

					if (tnext > tign[kl])
					{
						ts = tchar[k0];
						qdot[kl][0] = hb * std::max((tf - ts), 0.0);
					}
				}
			}
		}//Major l loop
	}//Major k loop

	//!c Update fractions ignited and burned out, to apply at next step start

	for (int k = 1; k <= numFuelTypes; k++)
	{
		int k0 = k - 1;

		flit[k0] = 0.0;
		fout[k0] = 0.0;
		
		for (int l = 0; l <= k; l++)
		{
			kl = Loc(k, l);
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

//Loc():
/** This function converts the indexes of the pairwise fuel interaction triangular matrix space
 * to indexes of the arrays used to represent it for several computed variables.
 *
 * @param[in] k Triangular matrix column (row) index, (1 - number of fuel types).
 * @param[in] l Triangular matrix row (column) index, (0 - k), = partner.
 *              This index starts at 0, which represent the "no companion" pairs.
 *
 * @returns The compact array index representing triangular matrix position [k, l], notated as kl. 
 * 
 * @par Indexes:
 * The matrix indices representing the fuel types and are 1 based, i.e. 1 = fuel type 1.  The two
 * dimensions of the matrix are represented by index variables k and l, where k represents the first
 * fuel type and the l index represents the partner fuel, with 0 representing no partner.  In this
 * way the indexes for the two dimensions have the same meaning for all positive values.  In the
 * orignal Fortran, which uses 1 based array indexing by default, the fuel type value align with all
 * arrays organized by fuel type.  However, things are a bit prone to confusion in C++.  There is no
 * easy way to make the triangular matrix more C-ish since one dimension is already 0 based.
 * Therefore we keep the indexing scheme the same.  Calling code must keep the index model in mind.
 * However, this function's return valued differes from the Fortran implementation in that we assume
 * the arrays used to represent triangular matrix data use native.
 *
 * The fact that fuel types start at 1 in triangular matrix space but start with 0 in C array space
 * can make the code confusing.  To make things clearer we have introduced the convention of using
 * k and l only where 1 based indexes are used.  Where the 0 based equivalent is needed k0, and to a
 * lesser extent l0????? are substituted.
 * 
 * @note: This will only return valid (occupied) coordinates of the triangular matrix.
 * The code assumes maxno and maxkl are the maximum dimensions.
 * Further error checking would require that the number of fuel classes be know. ?????
 *
 * @par History:
 * This function was originally implemented as a statement function defined in seven places in the
 * original Fortran code.
 *
 */
int Loc(const int k, const int l)
{
	int loc;//Return value: Index in a compact array representing the triangular matrix values (aka kl).

	//Input validity checking:
	if ((k < 1) || (k > NumFuelTypes))
	{
		Stop("Loc(): Invalid value of k " + std::to_string(k));
	}

	if ((l < 0) || (l > k))
	{
		Stop("Loc(): Invalid value of l " + std::to_string(l));
	}

	loc = k * (k + 1) / 2 + l;
	loc -= 1;//Convert to 0 based index.

	//Check value calculated is in the valid range:
	if ((loc < 0) || (loc > (Length_kl(NumFuelTypes) - 1)))
	{
		Stop("Loc(): Invalid index returned " + std::to_string(loc));
	}

	return loc;
}

//Lenkl():
/** Compute the triangular matrix size for a given number of fuel types.
 *
 * The original code set this as a fixed value notated maxkl.
 *
 * @param numFuels The number of fuel type elements.
 *
 * @returns The triangular matrix size.
 *
 * @par History:
 * Added for C++.  This calculation was performed several places in the original model.
 */
int Length_kl(int numFuels)
{
	//Add one to one dimension for the 'no companion' interaction element:
	int lenkl = numFuels * (numFuels + 1) / 2 + numFuels;
	return lenkl;
}

//ErrorApprox()
/** Approximate the error function.
 *
 * I believe this function is the approximation of the complementary error function as
 * described in Albini 1995 and Albini & Reinhardt 1995.  It was obtained from Hastings (et
 * al.) 1955, but I cann't identify the equation in that reference.
 *
 * @param[in] h
 * @param[in] theta Temperature rise required for the start of moisture loss.
 *
 * @par History:
 * This function was originally implemented as a statement function in DRYTIM(), as f().
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

//ValidateOutputVector():
/** This utility function checks that the output vector passed in is the correct size and resizes it
 * if necessary.
 *
 * Allow output vectors to be empty and resize them if needed.  If the correct size do nothing.
 * Treat other sizes a potential problem and warn about them.
 *
 * @param[in,out]	output		The output vector 
 * @param[in]		outputName	
 *
 * @returns Nothing.
 * 
 * @par History:
 * Added for C++ to reduce code repetition in Simulate().
 */
void ValidateOutputVector(std::vector<double>& output, const std::string outputName)
{
	int lenkl = Length_kl(NumFuelTypes);

	if (output.size() != lenkl)
	{
		if (!output.empty())
		{
			Warning(outputName + " has unexpected size " + std::to_string(output.size()) +
			        ". Erasing and resizing.");
		}

		output.assign(lenkl, 0.0);//The incoming values are ignored so resize() would be fine too.
	}
}

//SaveStateToFile():
/** Output the state of the simulation at the current timestep to file (if needed):
 * Sequential calls to this routine will produce a full history of the simulated fire.
 * Data is only saved when the SaveHistory setting is set to true.
 *
 * Successful output from this routine is treated as non-critical as it doesn't impact the
 * simulation process.  Output failures are reported but are not treated as fatal.
 *
 * The file name is currently fixed.  A history file from a previous run will prevent a new one
 * from being created.  In the future we may automatically number the file if it already exists.
 * Allowing the file name to be specified would only work in some contexts.  For example, we
 * currently have no way to pass in a file name via the R interface.
 *
 * @param[in] ts		Current timestep count.
 * @param[in] time		Current time (s).
 * @param[in] number	Actual number of fuel components.
 * @param[in] parts		Fuel component names / labels. [maxno]
 *
 * All the outputs from START() and STEP():
 * Currently only some of these are passed in to be saved. In the future others may be added.
 * @param[in] wo		Current ovendry loading for the larger of each component pair, kg / sq m. [maxkl]
 * @param[in] diam		Current diameter of the larger of each fuel component pair, m. [maxkl]	!!!!!
	!real*4, intent(in) :: flit(maxno)		Fraction of each component currently alight.
	!real*4, intent(in) :: fout(maxno)		Fraction of each component currently gone out.
 * The following are recomputed at each time point but only the final values would be needed:
	!real*4, intent(in) :: tdry(maxkl)		! Time of drying start of the larger of each fuel component pair, s.
	!real*4, intent(in) :: tign(maxkl)		! Ignition time for the larger of each fuel component pair, s. [maxkl]
	!real*4, intent(in) :: tout(maxkl)		! Burnout time of larger component of pairs, s.
	!real*4, intent(in) :: qcum(maxkl)		! Cumulative heat input to larger of pair, J / sq m.
	!real*4, intent(in) :: tcum(maxkl)		! Cumulative temp integral for qcum (drying).
	!real*4, intent(in) :: acum(maxkl)		! Heat pulse area for historical rate averaging.
 * This includes time so again only the final value would be needed:
	!real*4, intent(in) :: qdot(maxkl, mxstep)	! History (post ignite) of heat transfer rate
												! to the larger of each component pair, W / sq m..
	!real*4, intent(in) :: ddot(maxkl)		! Diameter reduction rate, larger of pair, m / s.
	!real*4, intent(in) :: wodot(maxkl)		! Dry loading loss rate for larger of pair.
	!real*4, intent(in) :: work(maxno)	! Workspace array.
 * @param[in] fi		Current fire intensity (site avg), kW / sq m.
 *
 * @returns Nothing.
 *
 * @par Format for the variable output:
 * The data is written in long format with tab delimited fields:
 * timestep (integer), time (float), variable (string), value (float), and IDs (strings).
 * The IDs are currently only used to identify the fuel pairs.  If that is the only use they
 * should be renamed.
 *
 * @par History:
 * Added for module.  This routine provides an alternative to the original STASH() function. [More...]
 */
void SaveStateToFile(const int ts, const double time, const int number,
                     const std::vector<std::string> parts, const std::vector<double> wo,
                     const std::vector<double> diam, const double fi)//const
{
	//Local constants:
	const std::string histFileName("BurnupHistory.txt");//histFile in Fortran.
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
		//for (int k0 = 0; k0 < number; k0++)
		//for (int k = 0; k <= number; k++)//Bad!!!!!
		for (int k = 1; k <= number; k++)
		{
			int k0 = k - 1;//Only used once.
			fuelName = parts[k0];

			//do l = 0, k
			//for (int l = 0; l < number; l++)
			for (int l = 0; l <= k; l++)
			{
				//int l0 = l - 1;
				//int kl = Loc(k0 + 1, l);//Triangular matrix index, 0 based.
				int kl = Loc(k, l);

				//Get the name of the partner component:
				if (l == 0)
				{
					compName = noCmpStr;
				}
				else
				{
					compName = parts[l - 1];//l0
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
