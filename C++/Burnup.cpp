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

#include <cmath>//<math.h>//For pow()/
#include <iostream>//Or just <ostream>?

#include "Burnup.h"

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



//InteractiveUI()

//Simulate()

//SimulateR()

//DUFBRN()

//DufBrnR()

//GETDAT()

//...

//RETRY()

//ARRAYS

//SORTER

//OVLAPS

//START

//FIRINT

//STASH

//TIGNIT

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
	int n;//Counter
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

//SUMMARY

//Loc
/**
! This function converts the indexes of the pairwise fuel interaction triangular matrix space
! to indexes of the arrays used to represent it for several computed variables.
!
! Note: This will only return valid (occupied) coordinates of the triangular matrix.
! Error checking would require that the number of fuel classes be know. 
!
! History: This function was originally implemented as a statement function defined in seven places
! in the original code.

 * @param k Triangular matrix column (row) index, (1:number of fuel types).
 * @param l Triangular matrix row (column) index, (0:k), = partner.
 *          This index starts at 0, which represent the "no companion" pairs.
 */
int Loc(const int k, const int l)
{
	int loc;//Return value: Index in a compact array representing the triangular matrix values.

	//Validity checking:
	if ((k < 1) || (k > maxno))
	{
		//print *, "Loc(): Invalid value of k ", k
		std::cout << "Loc(): Invalid value of k " << k << std::endl;
	}

	if ((l < 0) || (l > k))
	{
		//print *, "Loc(): Invalid value of l ", l
		std::cout << "Loc(): Invalid value of l " << l << std::endl;
	}

	loc = k * (k + 1) / 2 + l;

	if ((loc < 1) || (loc > maxkl))
	{
		//print *, "Loc(): Invalid index returned ", loc
		std::cout << "Loc(): Invalid index returned " << loc << std::endl;
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
