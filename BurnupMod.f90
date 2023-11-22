!---------------------------------------------------------------------------------------------------
! BurnupMod.f90
! Burnup Wildfire Fuel Consumption Module
!
! Original code by: Frank A. Albini
! Edited and updated by: Joshua M. Rady
! Woodwell Climate Research Center
! 9/2023
!
! This is an reimplementation of the Burnup wildfire fuel consumption model as a Fortran module.
! 
! 	This module is based on the Burnup model by Frank Albini and collaborators (see references).
! The original Burnup source code (Fortran IV/66/77, fixed form) was previously updated to modern
! Fortran (2003+).  Here that interactive executable has be reformulated into a form that can be
! compiled as a linkable module or shared library.
!
! 	The goal was to provide a version of the Burnup model that could be easily embedded in other
! code or coupled with other models.  The code produces identical results to the original code and
! has no dependancies.  The intent is that calling code only need to know Burnup's input and
! outputs.
!
! 	Routines providing C interoperable interfaces have been included for the main program entry
! points to allow the code to interface with R when compiled as a shared library (.so file).
!
! 	The routines that run the UI of the original program have been left for now but may be removed.
!
! Formating:
! 	I have used tabs for code indenting and for alignment of comments.  This is to aid with reading
! and future porting of this code (e.g. to C++).  While tabs are not valid Fortran characters,
! compilers tend to tolerate them.  I have used 4-space equivalent tabs for layout purposes here,
! and they can easily be converted in the future.
!
! In-code Documentation:
! 	Many comments have been added to the code to increase it's readability and to document changes.
! Comments from the original code are marked with '!c'.  Comments starting with just '!' have been
! added.
!
! 	Routines from the original propgram maintain their all caps format while new routines use
! camelcase.
!
! References:
! 	The original source code was obtained from
!
! Program BURNUP, a simulation model of the burning of large woody natural fuels.
! Frank A. Albini
! Missoula, MT: USDA Forest Service. Unpublished report. Research Grant INT-92754-GR. 1994.
! Appendix B.
!
! The original report is held by the National Forest Service Library in Fort Collins.
!
! ...
!
! Caveats:
! 	This code compiles, runs, and produces output identical to the original Burnup code in tests but
! the code is under active development and is not complete, so errors could be present that have not
! been identified.
!
! 	There is no licence provided for the original code.  It is though to be open by provenance, but
! that man not be correct.  The licence for this code is under consideration.
!---------------------------------------------------------------------------------------------------

!program BURNUP
module BurnupMod
	implicit none
	! private

	! These should be calculation only.
	public :: DUFBRN
	public :: ARRAYS
	public :: SORTER
	public :: OVLAPS
	public :: START
	public :: FIRINT
	public :: TIGNIT
	public :: DRYTIM
	public :: HEATX
	public :: TEMPF
	public :: STEP


	! Contain original menu based UI.  Probably not runnable in most contexts.
	! Should these just be removed?
	private :: BurnupMain
	private :: GETDAT
	private :: GetComponentParameters
	private :: GetFireAndSimProperties
	private :: ReviewDataMenu
	private :: AddDeleteComponent
	private :: ReviewFuelComp
	private :: ReviseFuelComponent
	private :: ReviseFuelParameters
	private :: ReviewFireEnvData
	private :: ReviewIntlCntlVars
	private :: GetDataFromFiles
	private :: ArchiveMenu
	private :: ArchiveSettings
	private :: RETRY
	
	! File IO:
	private :: STASH
	private :: SUMMARY
	
	! Utilities:
	private :: Loc
	private :: ErrorApprox
	private :: AskForReal1
	private :: AskForReal2
	private :: AskForReal

	! R wrappers:
	! These routines are intended to be only used from R not from Fortan.  As such they will only be
	! used when the code is compiled as a shared library rather than as part of an executable.
	! Specifying these as private still allows them to be used in the shared library context but
	! generates a compiler warning.
	public :: SimulateR
	!private :: DufBrnR
	public :: DufBrnR

	! Program level dimensional constants:
	! In the original code these parameters were passed into all routines that need.  Here we change
	! to module level scope.  This allows us to reduce the number of arguments to routines.
	! While currently fixed, these can probably be made dynamic to be set a initialization.
	
	! The maximum number of fuel components or types.  This is used to build fixed size data
	! structures.  The number of elements may be less than this so some indexes may be empty.
	! The original program fixes this arbitrarily at 10 fuel components.
	integer, parameter :: maxno = 12!10
	! The maximum number of non-zero entries in the triangular matrix of fuel interaction pairs:
	! Add one to one dimension for the  'no companion' interaction element.
	integer, parameter :: maxkl = maxno * (maxno + 1) / 2 + maxno
	! The maximum dimension of historical sequences (primarily for qdot):
	integer, parameter :: mxstep = 20

	! Physical constants:
	! In the original code these were defined as variables and shared through subroutine arguments.
	real*4, parameter :: ch2o = 4186.0	! Specific heat capacity of water, J / kg K
	real*4, parameter :: tpdry = 353.0	! Temperature (all components) start drying (K)

	!integer, parameter :: klVal = 7 ! JMR_TEMP_REPORTING: Loc(3,1) = 7
	!integer, parameter :: klVal = 11 ! JMR_TEMP_REPORTING: Loc(4,1) = 11
	!integer, parameter :: klVal = 29 ! JMR_TEMP_REPORTING: Loc(7,1) = 29
	integer, parameter :: klVal = -1 ! JMR_TEMP_REPORTING: This should disable reporting.

contains


	! The original contents of the main Burnup program.  At the moment this code has just beeen
	! moved here.  It can't be expected to work as is.
	subroutine BurnupMain()
		implicit none

		! Locals:
		! The original declarations in the original order:
		real*4 :: wdry(maxno)			! Ovendry mass loading, kg/sq m
		real*4 :: ash(maxno)			! Mineral content, fraction dry mass
		real*4 :: htval(maxno)			! Low heat of combustion, J / kg
		real*4 :: fmois(maxno)			! Moisture fraction of component
		real*4 :: dendry(maxno)			! Ovendry mass density, kg / cu m
		real*4 :: sigma(maxno)			! Surface to volume ratio, 1 / m
		real*4 :: cheat(maxno)			! Specific heat capacity, (J / K) / kg dry mass
		real*4 :: condry(maxno)			! Thermal conductivity, W / m K, ovendry
		real*4 :: alfa(maxno)			! Dry thermal diffusivity of component, sq m / s
		real*4 :: tpig(maxno)			! Ignition temperature, K
		real*4 :: tchar(maxno)			! Char temperature, K
		real*4 :: flit(maxno)			! Fraction of each component currently alight
		real*4 :: fout(maxno)			! Fraction of each component currently gone out
		real*4 :: work(maxno)			! Workspace array
		real*4 :: elam(maxno, maxno)	! Interaction matrix
		real*4 :: alone(maxno)			! Non-interacting fraction for each fuel class.
		real*4 :: area(maxno)			! Fraction of site area expected to be covered at
										! least once by initial planform area of ea size
		real*4 :: fint(maxno)			! Corrected local fire intensity for each fuel type.
		real*4 :: xmat(maxkl)			! Consolidated interaction matrix
		real*4 :: tdry(maxkl)			! Time of drying start of the larger of each
										! fuel component pair
		real*4 :: tign(maxkl)			! Ignition time for the larger of each
										! fuel component pair
		real*4 :: tout(maxkl) 			! Burnout time of larger component of pairs
		real*4 :: wo(maxkl)				! Initial dry loading by interaction pairs
		real*4 :: wodot(maxkl)			! Dry loading loss rate for larger of pair
		real*4 :: diam(maxkl)			! initial diameter, m [by interaction pairs]
		real*4 :: ddot(maxkl)  			! Diameter reduction rate, larger of pair, m / s
		real*4 :: qcum(maxkl) 			! Cumulative heat input to larger of pair, J / sq m
		real*4 :: tcum(maxkl) 			! Cumulative temp integral for qcum (drying)
		real*4 :: acum(maxkl) 			! Heat pulse area for historical rate averaging
		real*4 :: qdot(maxkl, mxstep)	! History (post ignite) of heat transfer rate
		integer :: key(maxno)			! Ordered index list
		character*12 :: parts(maxno)	! Fuel component names / labels
		character*12 :: list(maxno)		!
		character*12 :: infile			! Stores the name of input data files.
		character*12 :: outfil			! The name of the summary file.  This currently static in the code
										! below.  It would be better to ask.
		logical :: nohist				! Flag indicating if history output should be not be stored.

		! The rest in order of appearance:
		real :: fimin			! Fire intensity (kW / sq m) at which fire goes out. (Change to parameter?)
		integer :: nruns		! The number of simulations run.
		integer :: number		! The actual number of fuel classes
		real*4 :: fi			! Site avg fire intensity (kW / sq m)
		real*4 :: ti 			! Spreading fire residence time (s)
		real*4 :: u				! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4 :: d				! Fuelbed depth (m)
		real*4 :: tpamb			! Ambient temperature (K)
		real*4 :: ak			! Area influence factor (ak / K_a parameter)
		real*4 :: r0			! Minimum value of mixing parameter
		real*4 :: dr			! Max - min value of mixing parameter
		real*4 :: dt			! Time step for integration of burning rates (s)
		integer :: ntimes		! Number of time steps.
		real*4 :: wdf			! Duff loading (kg/m^2, aka W sub d)
		real*4 :: dfm			! Ratio of moisture mass to dry organic mass /
								! duff fractional moisture (aka R sub M).
		integer :: ihist		! User input value.
		real :: dfi 			! Duff fire intensity (aka I sub d) for DUFBRN().
		real :: tdf 			! Burning duration (aka t sub d) for DUFBRN().
		integer :: now			! Index marking the current time step.
		real :: tis				! Current time (ti + number of time steps * dt).
		!real*4 :: tpdry			! Temperature (all components) start drying (K)
		integer :: ncalls		! Counter of calls to this START().
		!real*4 :: ch2o			! Specific heat capacity of water, J / kg K
		real :: fid				! Fire intensity due to duff burning.
		integer :: nun			! Stash file unit identifier.
		integer :: readStat		! IO error status.
		integer :: next 		! User menu selection.

		! Constants:
		! STASH() leaves the history file open when it returns.  The original code hard codes this file
		! unit value.  It would be more robust if it were returned like nun.
		integer, parameter :: mum = 66 ! History file unit number.

		fimin = 0.1
		nruns = 0

		do
			call GETDAT(infile, outfil, parts, wdry, ash, htval, &
						fmois, dendry, sigma, cheat, condry, tpig, &
						tchar, number, maxno, fi, ti, u, d, tpamb, &
						ak, r0, dr, dt, ntimes, wdf, dfm, nruns, &
						area, fint)

			do
				write(*, "(' Enter 1 to store fire history, 0 to skip ' ,$)")
				read(*, *) ihist
				! Could add read error checking here.
		
				! If a valid selection was made continue, otherwise prompt again:
				if ((ihist .eq. 0) .or. (ihist .eq. 1)) then
					nohist = ihist .eq. 0
					exit
				end if
			end do

			call ARRAYS(maxno, number, wdry, ash, dendry, fmois, &
						sigma, htval, cheat, condry, tpig, tchar, &
						diam, key, work, ak, elam, alone, xmat, wo, &
						maxkl, parts, list, area)
			! Note: elam and alone are passed in and returned but are not used again.

			call DUFBRN(wdf, dfm, dfi, tdf)

			! Initialize variables and data structures:
			now = 1
			tis = ti
			call START(tis, mxstep, now, maxno, number, wo, alfa, &
						dendry, fmois, cheat, condry, diam, tpig, &
						tchar, xmat, tpamb, fi, flit, fout, &
						tdry, tign, tout, qcum, tcum, acum, qdot, &
						ddot, wodot, work, u, d, r0, dr, &
						ncalls, maxkl)

			! If the duff burns longer than the passing fire front then have it's intensity
			! contribute to the post front fire environment, otherwise ignore it:
			if (tis .lt. tdf) then
				fid = dfi
			else
				fid = 0.0
			end if

			if (nohist .eqv. .false.) then
				call STASH(tis, now, maxno, number, outfil, fi, flit, &
							fout, wo, wodot, diam, ddot, tdry, tign, &
							tout, fmois, maxkl, nun)
			end if

			! Calculate the initial fire intensity:
			call FIRINT(wodot, ash, htval, maxno, number, maxkl, area, fint, fi)

			! If the fire intensity is above the extinguishing threshold calculate combustion until
			! the fire goes out or the number of timesteps is reached:
			if (fi .gt. fimin) then
				do while (now .LT. ntimes)
					call STEP(dt, mxstep, now, maxno, number, wo, alfa, &
								dendry, fmois, cheat, condry, diam, tpig, &
								tchar, xmat, tpamb, fi, flit, fout, &
								tdry, tign, tout, qcum, tcum, acum, qdot, &
								ddot, wodot, work, u, d, r0, dr, &
								ncalls, maxkl, tis, fint, fid)
	
					! Update time trackers:
					now = now + 1
					tis = tis + dt
					
					! Get the duff contribution while it remains burning:
					if (tis .lt. tdf) then
						fid = dfi
					else
						fid = 0.0
					end if
			
					! Calculate the fire intensity at this time step:
					call FIRINT(wodot, ash, htval, maxno, number, maxkl, area, fint, fi)

					if (fi .LE. fimin) then
						exit
					else
						if (nohist .eqv. .false.) then
							call STASH(tis, now, maxno, number, outfil, fi, flit, &
										fout, wo, wodot, diam, ddot, tdry, tign, &
										tout, fmois, maxkl, nun)
						end if
					end if
				end do
			end if ! (fi .gt. fimin)

			close(nun)

			outfil = 'SUMMARY.DAT'
			call SUMMARY(outfil, number, maxno, maxkl, parts, nun, &
							tis, ak, wdry, fmois, sigma, tign, tout, xmat, wo, diam)
			close(nun)

			do
				write(*, "(' Exercise completed. Do another = 1 , quit = 0 ',$)")
				read(*, *, iostat = readStat) next
				if (readStat .eq. 0) then

					if (next .eq. 0) then
						close(mum)
						stop' Terminated'
					else if (next .eq. 1) then
						exit
					!Otherwise an invalid value was provided. Ask again.
					end if
				end if ! (readStat .eq. 0)
			end do

			close(mum)
		end do ! Main program loop.

	end subroutine BurnupMain


	! Perform a simulation with prescribed inputs and return fuel consumption properties.
	!
	! History: Added as an programatic alternative entry point to the original interactive program.
	!
	! ToDo:
	! - Provide a way to specify if fire history should be stored and returned.
	subroutine Simulate(fi, ti, u, d, tpamb, ak, r0, dr, dt, wdf, dfm, ntimes, number, &
						parts, wdry, ash, htval, fmois, dendry, sigma, cheat, condry, tpig, tchar, &
						xmat, tign, tout, wo, diam)
		implicit none

		! Arguments:
		! Igniting fire and environmental data:
		
		! The value passed in for fi is the fire front intensity.  The variable is later reused and
		! updated by FIRINT().  It is passed on to other routines that use bu do not change it.
		! These two uses could be separated.  The value returned the value, is the final intensity,
		! which might be of use.  A history would be more valuable.
		real*4, intent(inout) :: fi		! Current fire intensity (site avg), kW / sq m
		real*4, intent(in) :: ti		! Igniting fire residence time (s).
		real*4, intent(in) :: u			! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4, intent(in) :: d			! Fuelbed depth (m)
		real*4, intent(in) :: tpamb		! Ambient temperature (K)
		! Internal and control variables:
		real*4, intent(in) :: ak		! Area influence factor (ak / K_a parameter)
		real*4, intent(in) :: r0		! Minimum value of mixing parameter
		real*4, intent(in) :: dr		! Max - min value of mixing parameter
		!real*4, intent(in) :: dt		! Time step for integration of burning rates (s)
		real*4, intent(inout) :: dt		! Time step for integration of burning rates (s).
										! On completion contains the time the fire went out.
		
		! Considering removing these two.  See below:
		real*4, intent(in) :: wdf		! Duff loading (kg/m^2, aka W sub d)
		real*4, intent(in) :: dfm		! Ratio of moisture mass to dry organic mass /
										! duff fractional moisture (aka R sub M).
		integer, intent(in) :: ntimes	! Number of time steps.  Move down?
		integer, intent(in) :: number	! The number of fuel classes. ! Could try to remove?????
		
		! Fuel component property arrays:  The values will not change but they may be reordered.
		! Returning the reordered arrays may be overkill.  The revised order might be sufficient.
		! However, setting these to inout allows the values to be modified for use interally in this
		! routine, which is useful for now.
		! JMR_NOTE: Some of these, e.g. SAV can safely be made in only!!!!!
		
		! Character strings can't be passed in from R so we leave part blank for now (see below):
		character*12, intent(inout) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(inout) :: wdry(maxno)		! Ovendry mass loading, kg/sq m
		real*4, intent(inout) :: ash(maxno)			! Mineral content, fraction dry mass
		real*4, intent(inout) :: htval(maxno)		! Low heat of combustion, J / kg,					AKA heat content
		real*4, intent(inout) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(inout) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(inout) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(inout) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass		J / kg K
		real*4, intent(inout) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		
		
		real*4, intent(inout) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(inout) :: tchar(maxno)		! Char temperature, K
		
		!integer, intent(in) :: maxno				! The maximum number of fuel classes allowed.
		! Could try to remove?:
		!integer, intent(in) :: number			! The number of fuel classes.  Move up?
		
		! Calculated outputs:
		! The following are the main variables output by SUMMARY(): [name], fr, ti, to, wd, di
		real*4, intent(out) :: xmat(maxkl)			! Table of influence fractions between components
		real*4, intent(out) :: tign(maxkl)			! Ignition time for the larger of each fuel component pair
		real*4, intent(out) :: tout(maxkl)			! Burnout time of larger component of pairs
		real*4, intent(out) :: wo(maxkl)			! Current ovendry loading for the larger of
													! each component pair, kg/sq m
		real*4, intent(out) :: diam(maxkl)			! Current diameter of the larger of each
													! fuel component pair, m
		
		! Locals:
		! Arrays:
		real*4 :: alfa(maxno)			! Dry thermal diffusivity of component, sq m / s
		real*4 :: flit(maxno)			! Fraction of each component currently alight
		real*4 :: fout(maxno)			! Fraction of each component currently gone out
		real*4 :: work(maxno)			! Workspace array
		real*4 :: elam(maxno, maxno)	! Interaction matrix
		real*4 :: alone(maxno)			! Non-interacting fraction for each fuel class.
		real*4 :: area(maxno)			! Fraction of site area expected to be covered at
										! least once by initial planform area of ea size
		real*4 :: fint(maxno)			! Corrected local fire intensity for each fuel type.
		real*4 :: tdry(maxkl)			! Time of drying start of the larger of each
		real*4 :: wodot(maxkl)			! Dry loading loss rate for larger of pair
		real*4 :: ddot(maxkl)  			! Diameter reduction rate, larger of pair, m / s
		real*4 :: qcum(maxkl) 			! Cumulative heat input to larger of pair, J / sq m
		real*4 :: tcum(maxkl) 			! Cumulative temp integral for qcum (drying)
		real*4 :: acum(maxkl) 			! Heat pulse area for historical rate averaging
		real*4 :: qdot(maxkl, mxstep)	! History (post ignite) of heat transfer rate
		integer :: key(maxno)			! Ordered index list
		! Temporary: Leave the fuel component names blank.
		!character*12 :: parts(maxno)	! Fuel component names / labels
		character*12 :: list(maxno)		!
		!logical :: nohist				! Flag indicating if history output should be not be stored.
		
		! Scalars in order of appearance:
		real :: dfi 			! Duff fire intensity (aka I sub d) for DUFBRN().
		real :: tdf 			! Burning duration (aka t sub d) for DUFBRN().
		integer :: now			! Index marking the current time step.
		real :: tis				! Current time (ti + number of time steps * dt).
		integer :: ncalls		! Counter of calls to START().
		real :: fid				! Fire intensity due to duff burning.
		
		
		! In the original code fmin was a local treated as a constant.  Passing it in might be good:
		real, parameter :: fimin = 0.1 ! Fire intensity (kW / sq m) at which fire goes out.

		! There are a large number of locals in this routine that are not explictly initialized.
		! Most are initialized in ARRAYS() and START.  Testing was done to confrim that explicit
		! initialization was not needed for the remaining locals.

		! JMR: Temporary:
		! Parts will be empty:
! 		parts = "Test"
! 		call ArchiveSettings(parts, wdry, ash, htval, fmois, dendry, &
! 								sigma, cheat, condry, tpig, tchar, maxno, number, &
! 								fi, ti, u, d, tpamb, &
! 								ak, r0, dr, dt, wdf, dfm, ntimes)
! 		
! 		return
		
		
		!print *, "ARRAYS:" ! JMR: Temp reporting.
		! Sort the fuel components and calculate the interaction matrix...
		call ARRAYS(maxno, number, wdry, ash, dendry, fmois, &
					sigma, htval, cheat, condry, tpig, tchar, &
					diam, key, work, ak, elam, alone, xmat, wo, &
					maxkl, parts, list, area)

		!print *, "DUFBRN:" ! JMR: Temp reporting.
		! The original code calls DUFBRN() here.  I'm leaving this here while getting the code
		! running but it would probably better to pass the output of DUFBRN() in instead.
		call DUFBRN(wdf, dfm, dfi, tdf)
		
		! Initialize variables and data structures:
		! Which of these variables need to be declared?????
		now = 1
		tis = ti
		!print *, "START:" ! JMR: Temp reporting.
		call START(tis, mxstep, now, maxno, number, wo, alfa, &
					dendry, fmois, cheat, condry, diam, tpig, &
					tchar, xmat, tpamb, fi, flit, fout, &
					tdry, tign, tout, qcum, tcum, acum, qdot, &
					ddot, wodot, work, u, d, r0, dr, &
					ncalls, maxkl)
		
		! If the duff burns longer than the passing fire front then have it's intensity
		! contribute to the post front fire environment, otherwise ignore it:
		if (tis .lt. tdf) then
			fid = dfi
		else
			fid = 0.0
		end if
		
		! Record the state at this timepoint somehow...
		
		
		! Calculate the initial fire intensity:
		!print *, "FIRINT:" ! JMR: Temp reporting.
		call FIRINT(wodot, ash, htval, maxno, number, maxkl, area, fint, fi)
		
		!print *, "Loop:" ! JMR: Temp reporting.
		! If the fire intensity is above the extinguishing threshold calculate combustion until
		! the fire goes out or the number of timesteps is reached:
		if (fi .gt. fimin) then
			do while (now .LT. ntimes)
				!print *, "STEP:" ! JMR: Temp reporting.
				call STEP(dt, mxstep, now, maxno, number, wo, alfa, &
							dendry, fmois, cheat, condry, diam, tpig, &
							tchar, xmat, tpamb, fi, flit, fout, &
							tdry, tign, tout, qcum, tcum, acum, qdot, &
							ddot, wodot, work, u, d, r0, dr, &
							ncalls, maxkl, tis, fint, fid)

				
				! JMR: Temp reporting:
				if (now .eq. 1 .or. mod(now, 100) .eq. 0) then
					!print *, "Next"
				end if
				
				! Update time trackers:
				now = now + 1
				tis = tis + dt
				
				! Get the duff contribution while it remains burning:
				if (tis .lt. tdf) then
					fid = dfi
				else
					fid = 0.0
				end if
		
				! Calculate the fire intensity at this time step:
				call FIRINT(wodot, ash, htval, maxno, number, maxkl, area, fint, fi)

				if (fi .LE. fimin) then
					exit
				else
					! Record the state at this timepoint somehow...
				end if
			end do
		end if ! (fi .gt. fimin)

		! The main fire properties are returned as output arguments...
		! A replacement for SUMMARY() could go here.
		! Return the time when the fire dropped below fimin as the time the fire "went out".
		! We return the value in dt.  This may be changed in the future.
		dt = tis
		! We could also return the number of timesteps completed in ntimes but that doesn't add much.

	end subroutine Simulate


	! This is a wrapper for Simulate() that allows it to be called from R:
	! 
	! History: Added for module.
	subroutine SimulateR(fi, ti, u, d, tpamb, ak, r0, dr, dt, wdf, dfm, ntimes, number, &
							wdry, ash, htval, fmois, dendry, sigma, cheat, condry, tpig, tchar, &
							xmat, tign, tout, wo, diam) bind(C, name = "simulater")
		implicit none

		! Arguments:
		double precision, intent(inout) :: fi		! Current fire intensity (site avg), kW / sq m
		double precision, intent(in) :: ti			! Igniting fire residence time (s).
		double precision, intent(in) :: u			! Mean horizontal windspeed at top of fuelbed (m/s).
		double precision, intent(in) :: d			! Fuelbed depth (m)
		double precision, intent(in) :: tpamb		! Ambient temperature (K)
		double precision, intent(in) :: ak			! Area influence factor (ak / K_a parameter)
		double precision, intent(in) :: r0			! Minimum value of mixing parameter
		double precision, intent(in) :: dr			! Max - min value of mixing parameter
		double precision, intent(inout) :: dt		! Time step for integration of burning rates (s)
		double precision, intent(in) :: wdf			! Duff loading (kg/m^2, aka W sub d)
		double precision, intent(in) :: dfm			! Ratio of moisture mass to dry organic mass /
													! duff fractional moisture (aka R sub M).
		integer, intent(in) :: ntimes						! Number of time steps.
		integer, intent(in) :: number	! The number of fuel classes. ! Could try to remove?????
		
		double precision, intent(inout) :: wdry(maxno)		! Ovendry mass loading, kg/sq m
		double precision, intent(inout) :: ash(maxno)		! Mineral content, fraction dry mass
		double precision, intent(inout) :: htval(maxno)		! Low heat of combustion, J / kg
		double precision, intent(inout) :: fmois(maxno)		! Moisture fraction of component
		double precision, intent(inout) :: dendry(maxno)	! Ovendry mass density, kg / cu m
		double precision, intent(inout) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		double precision, intent(inout) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		double precision, intent(inout) :: condry(maxno)	! Thermal conductivity, W / m K, ovendry
		double precision, intent(inout) :: tpig(maxno)		! Ignition temperature, K
		double precision, intent(inout) :: tchar(maxno)		! Char temperature, K
		
		! Calculated outputs:
		double precision, intent(out) :: xmat(maxkl)		! Table of influence fractions between components
		double precision, intent(out) :: tign(maxkl)		! Ignition time for the larger of each fuel component pair
		double precision, intent(out) :: tout(maxkl)		! Burnout time of larger component of pairs
		double precision, intent(out) :: wo(maxkl)			! Current ovendry loading for the larger of
															! each component pair, kg/sq m
		double precision, intent(out) :: diam(maxkl)		! Current diameter of the larger of each
															! fuel component pair, m

		! Local type conversion intermediates:
		real :: fiReal, dtReal
		real, dimension(maxno) :: wdryOut, ashOut, htvalOut, fmoisOut, dendryOut, sigmaOut
		real, dimension(maxno) :: cheatOut, condryOut, tpigOut, tcharOut
		real, dimension(maxkl) :: xmatR, tignR, toutR, woR, diamR

		integer :: i ! Counter
		character*12 :: parts(maxno) = ""	! Fuel component names / labels

		!print *, "tout", tout ! JMR_TEMP_REPORTING
		!print *, "toutR", toutR ! JMR_TEMP_REPORTING
		do i = 1, 3
			!parts(i) = "Fuel " // char(i)
			print *, char(i)
			print *, "Fuel " // char(i)
		end do

		fiReal = real(fi)
		dtReal = real(dt)
		wdryOut = real(wdry)
		ashOut = real(ash)
		htvalOut = real(htval)
		fmoisOut = real(fmois)
		dendryOut = real(dendry)
		sigmaOut = real(sigma)
		cheatOut = real(cheat)
		condryOut = real(condry)
		tpigOut = real(tpig)
		tcharOut = real(tchar)

		call Simulate(fiReal, real(ti), real(u), real(d), real(tpamb), &
						!real(ak), real(r0), real(dr), real(dt), &
						real(ak), real(r0), real(dr), dtReal, &
						real(wdf), real(dfm), ntimes, number, &
						parts, &
						wdryOut, ashOut, htvalOut, fmoisOut, dendryOut, &
						sigmaOut, cheatOut, condryOut, tpigOut, tcharOut, &
						xmatR, tignR, toutR, woR, diamR)

		fi = dble(fiReal)
		dt = dble(dtReal)
		wdry = dble(wdryOut)
		ash = dble(ashOut)
		htval = dble(htvalOut)
		fmois = dble(fmoisOut)
		dendry = dble(dendryOut)
		sigma = dble(sigmaOut)
		cheat = dble(cheatOut)
		condry = dble(condryOut)
		tpig = dble(tpigOut)
		tchar = dble(tcharOut)
		
		xmat = dble(xmatR)
		tign = dble(tignR)
		tout = dble(toutR)
		wo = dble(woR)
		diam = dble(diamR)
		
		!print *, "toutR", toutR ! JMR_TEMP_REPORTING
		!print *, "tout", tout ! JMR_TEMP_REPORTING

	end subroutine SimulateR


	!c Duff burning rate (ergo, intensity) and duration
	!
	! History: Modernized original Burnup subroutine.
	subroutine DUFBRN(wdf, dfm, dfi, tdf)
		implicit none

		! Arguments:
		real, intent(in) :: wdf		! Duff loading (kg/m^2, aka W sub d)
		real, intent(in) :: dfm 	! Ratio of moisture mass to dry organic mass /
									! duff fractional moisture (aka R sub M)
		real, intent(out) :: dfi	! Duff fire intensity (aka I sub d)
		real, intent(out) :: tdf	! Burning duration (aka t sub d)

		! Locals:
		real ff ! Fractional duff reduction depth from Brown et al. 1985, (aka F in report equation 4)

		dfi = 0.0
		tdf = 0.0
		if ((wdf .LE. 0.) .OR. (dfm .GE. 1.96)) return
		dfi = 11.25 - 4.05 * dfm
		ff = 0.837 - 0.426 * dfm
		tdf = 1.e+04 * ff * wdf / (7.5 - 2.7 * dfm)
	end subroutine DUFBRN


	! This is a wrapper for DUFBRN() that allows it to be called from R:
	! 
	! History: Added for module.
	subroutine DufBrnR(wdf, dfm, dfi, tdf) bind(C, name = "dufbrnc")
		implicit none

		! Arguments:
		double precision, intent(in) :: wdf		! Duff loading (kg/m^2, aka W sub d)
		double precision, intent(in) :: dfm 	! Ratio of moisture mass to dry organic mass /
												! duff fractional moisture (aka R sub M)
		double precision, intent(out) :: dfi	! Duff fire intensity (aka I sub d)
		double precision, intent(out) :: tdf	! Burning duration (aka t sub d)

		! Local type conversion intermediates:
		real :: dfiOut
		real :: tdfOut

		call DUFBRN(real(wdf), real(dfm), dfiOut, tdfOut)

		dfi = dble(dfiOut)
		tdf = dble(tdfOut)

	end subroutine DufBrnR


! -- Pagebreak --
! Pg. 78:


	!c This routine prompts user for input data in the form of file names
	!c or the direct input of numerical quantities.
	!c Routine makes consistency checks, tests nominal quantities, and
	!c allows review / revision of data before starting program,
	!c and allows archiving of run conditions for future use.
	!
	! History:
	! The original GETDAT() code used a lot of convoluted goto logic that was hard to follow.
	! This code has been refactored to a set of nested routines for improved readability.
	! Testing status: Menu reproduces original behavior.
	subroutine GETDAT(infile, outfil, parts, wdry, ash, htval, &
						fmois, dendry, sigma, cheat, condry, tpig, &
						tchar, number, maxno, fi, ti, u, d, tpamb, &
						ak, r0, dr, dt, ntimes, wdf, dfm, nruns, &
						area, fint)
		implicit none

		! Arguments:
		! JMR_NOTE: All the data variables are intent(out), which I thinks is right, even on a second call.
		! However, code below is sugests inout might be right.
		character*12, intent(out) :: infile			! Stores the name of input data files.
		! NOTE: There seems to be no reason that infile is returned.  It doesn't appear to be used,
		! and could be safely omitted.  It is maintained only to preserve the original interface.
		! This updated version of GETDAT() does not update infile. To reproduce the original
		! behavior and return the value outfil would need to be passed back by GetDataFromFiles().
		character*12, intent(out) :: outfil			! The name of the output file.  [Not used.]
		! Note: outfil is not actually used here and is maintained only to preserve the original
		! interface.  In the original code this argument is not initialized before being being
		! passed in, which clearly implies it is to be returned.  However, examination of the code
		! reveals that while the value is subsequently passed to other routines, its value is always
		! overridden.  Therefore it could be safely omitted.  This updated version of GETDAT() does
		! not update outfil.  To reproduce the original behavior and return the value outfil would
		! need to be passed back down the following call chain:
		! ArchiveMenu() -> ArchiveMenu() -> ReviewDataMenu -> GETDAT()
		character*12, intent(out) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(out) :: wdry(maxno) 			! Ovendry mass loading, kg/sq m
		real*4, intent(out) :: ash(maxno)			! Mineral content, fraction dry mass
		real*4, intent(out) :: htval(maxno)			! Low heat of combustion, J / kg
		real*4, intent(out) :: fmois(maxno)			! Moisture fraction of component
		real*4, intent(out) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(out) :: sigma(maxno)			! Surface to volume ratio, 1 / m
		real*4, intent(out) :: cheat(maxno)			! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(out) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(out) :: tpig(maxno)			! Ignition temperature, K
		real*4, intent(out) :: tchar(maxno)			! Char temperature, K
		integer, intent(out) :: number				! The actual number of fuel classes.
		integer, intent(in) :: maxno				! The maximum number of fuel classes allowed.
		real*4, intent(out) :: fi					! Current fire intensity (site avg), kW / sq m
		real*4, intent(out) :: ti					! Igniting fire residence time (s)
		real*4, intent(out) :: u					! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4, intent(out) :: d					! Fuelbed depth (m)
		real*4, intent(out) :: tpamb				! Ambient temperature (K)
		real*4, intent(out) :: ak					! Area influence factor (ak / K_a parameter)
		real*4, intent(out) :: r0					! Minimum value of mixing parameter
		real*4, intent(out) :: dr					! Max - min value of mixing parameter
		real*4, intent(out) :: dt					! Time step for integration of burning rates (s)
		integer, intent(out) :: ntimes				! Number of time steps.
		real*4, intent(out) :: wdf					! Duff loading (kg/m^2, aka W sub d)
		real*4, intent(out) :: dfm					! Ratio of moisture mass to dry organic mass /
													! duff fractional moisture (aka R sub M).
		integer, intent(out) :: nruns				! The number of times this routine as been run,
													! which controls variable initialization.
													! JMR_NOTE: There is currently no real reason
													! for this argument.  It is only used if more
													! than one run is performed in an interactive
													! session.  However, it doesn't need to be
													! passed back for this.  nruns could be made a
													! persistant local variable instead.
		real*4, intent(out) :: area(maxno)			! Fraction of site area expected to be covered at
													! least once by initial planform area of ea size
		real*4, intent(out) :: fint(maxno)			! Corrected local fire intensity for each fuel type.

		! Locals:
		integer :: i ! Counter
		integer :: in ! User menu selection.
		integer :: nn ! Component counter
		integer :: readStat ! IO error status.
		integer :: numb ! User menu selection.
		logical, save :: firstCall = .true. ! Is this the first time GETDAT() has been called?

		! If this is the first call to this routine, initialize the data:
		! JME_NOTE: This implies they should all be inout?
		if (nruns .eq. 0) then
			do i = 1 , maxno
				wdry(i) = 0.0
				ash(i) = 0.0
				htval(i) = 0.0
				fmois(i) = 0.0
				dendry(i) = 0.0
				sigma(i) = 0.0
				cheat(i) = 0.0
				condry(i) = 0.0
				tpig(i) = 0.0
				tchar(i) = 0.0
				area(i) = 0.0
				fint(i) = 0.0
				parts(i) = 'No data for'
			end do
			number = 0


! -- Pagebreak --
! Pg. 79:


			fi = 0.0
			ti = 0.0
			u = 0.0
			d = 0.0
			tpamb = 0.0
			ak = 0.0
			r0 = 0.0
			dr = 0.0
			dt = 0.0
			ntimes = 0
			wdf = 0.0
			dfm = 2.0
		end if

		! Determine how to obtain data:
		do ! Could convert to while with readStat = -1 at start.
			write(*, *) 'Enter 1 for keyboard entry of data'
			write(*, *) 'Enter 2 to fetch data stored in files'
			write(*, "(' Enter 0 to terminate program now ' , $)")

			read(*, *, iostat = readStat) in

			if (readStat .eq. 0) then
				if ((in .lt. 0) .or. (in .gt. 2)) then ! Invalid menu value.
					cycle
				else
					exit
				end if
			end if
		end do

		! Valid menu selections:
		if (in .eq. 0) then
			stop ' Program ended'
		else if (in .EQ. 2) then
			call GetDataFromFiles(parts, wdry, ash, htval, &
									fmois, dendry, sigma, cheat, condry, tpig, &
									tchar, number, maxno, &
									fi, ti, u, d, tpamb, &
									ak, r0, dr, dt, ntimes, wdf, dfm)
		else if (in .eq. 1) then ! User enters data...

			nruns = nruns + 1

			! Get the number of fuel elements:
			do
				write(*, "(' How many [ no more than'i3' ] fuel components ?    ',$)") maxno 
				read(*, *, iostat = readStat) numb
				if (readStat .eq. 0) then 
					if ((numb .LE. 0) .OR. (numb .GT. maxno)) then
						cycle
					else
						exit
					end if
				end if
			end do
			number = numb

			! For each element get the fuel parameters:
			nn = 0
			do while (nn .LT. number)
				nn = nn + 1

				call GetComponentParameters(parts, &
									wdry, ash, htval, fmois, dendry, &
									sigma, cheat, condry, tpig, tchar, nn)
			end do ! Fuel parameters.


! -- Pagebreak --
!  Pg. 80: All refactored.
! -- Pagebreak --
! Pg. 81:


			! Prompt for fire properties and simulation properties, but only the first time:
			if (firstCall) then
				call GetFireAndSimProperties(fi, ti, u, d, tpamb, ak, &
												r0, dr, dt, ntimes, wdf, &
												dfm)
				!firstCall = .false.
			end if

		end if ! (in .eq. X)

! -- Pagebreak --
! Pg. 82:

		!c All data collected --- now time to review and revise as needed
		call ReviewDataMenu(parts, wdry, ash, htval, fmois, dendry, &
							sigma, cheat, condry, tpig, tchar, maxno, number, &
							fi, ti, u, d, tpamb, &
							ak, r0, dr, dt, wdf, dfm, ntimes)


! -- Pagebreak --
! Pg. 83: All refactored.
! -- Pagebreak --
! Pg. 84: All refactored.
! -- Pagebreak --
! Pg. 85: All refactored.
! -- Pagebreak --
! Pg. 86: All refactored.
! -- Pagebreak --
! Pg. 87:


		firstCall = .false.

	end subroutine GETDAT


	! Request fuel parameters from the user for the specified fuel component (index = nFuel) and
	! store them in the passed data arrays.  Only this data index is stored.  The other indexes
	! and number of elements is not altered here.
	!
	! History: This routine reproduces code originally in GETDAT() (block starting at label 05).
	! Note: It would be nice to store the data arguments in an object.
	! Testing status: Initial tests reproduce original behavior.
	subroutine GetComponentParameters(parts, wdry, ash, htval, fmois, dendry, &
										sigma, cheat, condry, tpig, tchar, nFuel)
		implicit none

		! Arguments:
		character*12, intent(inout) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(inout) :: wdry(maxno)		! Ovendry mass loading, kg/sq m
		real*4, intent(inout) :: ash(maxno)			! Mineral content, fraction dry mass
		real*4, intent(inout) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(inout) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(inout) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(inout) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(inout) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(inout) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(inout) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(inout) :: tchar(maxno)		! Char temperature, K
		integer, intent(in) :: nFuel				! The fuel component to get (nn in original code).

		! Locals:
		real :: tigi, tchi ! Temperate intermediates in Celcius.
		character*12 :: name1

		! Constants:
		real, parameter :: small = 1.e-06
		real, parameter :: big = 1.e+06
		real, parameter :: ash1 = 0.0001
		real, parameter :: ash2 = 0.1
		real, parameter :: htv1 = 1.0e+07
		real, parameter :: htv2 = 3.0e+07
		real, parameter :: fms1 = 0.01
		real, parameter :: fms2 = 3.0
		real, parameter :: den1 = 200.0
		real, parameter :: den2 = 1000.0
		real, parameter :: sig1 = 4.0
		real, parameter :: sig2 = 1.0e+04
		real, parameter :: cht1 = 1000.0
		real, parameter :: cht2 = 2000.0 
		real, parameter :: con1 = 0.025
		real, parameter :: con2 = 0.25
		real, parameter :: tig1 = 200.0
		real, parameter :: tig2 = 400.0
		real, parameter :: tch1 = 250.0
		real, parameter :: tch2 = 500.0

		write(*, "(//' Fuel component number'i3 ' properties list')") nFuel
		write(*, "(' Name [ 12 characters or fewer ]    ', $)")
		read(*, "(a12)") name1
		parts(nFuel) = name1

		wdry(nFuel) = AskForReal1(name1, 'dry loading, kg / sq m ', small, big)

		ash(nFuel) = AskForReal1(name1, 'total mineral ash fraction    ', ash1, ash2)

		htval(nFuel) = AskForReal1(name1, 'heat of combustion, J / kg    ', htv1, htv2)

		fmois(nFuel) = AskForReal1(name1, 'moisture content, fraction    ', fms1, fms2)

		dendry(nFuel) = AskForReal1(name1, 'ovendry mass density, kg / cubic m    ', den1, den2)

		sigma(nFuel) = AskForReal1(name1, 'surface / volume ratio , inverse m    ', sig1, sig2)

		cheat(nFuel) = AskForReal1(name1, 'specific heat capacity ovendry , J / kg deg K  ', cht1, cht2)

		condry(nFuel) = AskForReal1(name1, 'ovendry thermal conductivity , W / m deg K    ', con1, con2)

		tigi = AskForReal1(name1, 'ignition temperature, deg C    ', tig1, tig2)
		tpig(nFuel) = tigi + 273.0

		tchi  = AskForReal1(name1, 'char [end pyrolysis] temperature, deg C    ', tch1, tch2)
		tchar(nFuel) = tchi + 273.0

	end subroutine GetComponentParameters


	! Get the fire and simulation parameters from the user:
	!
	! History: This routine reproduces code originally in GETDAT() (block starting at label 29).
	! Testing status: Initial tests reproduce original behavior.
	subroutine GetFireAndSimProperties(fi, ti, u, d, tpamb, ak, &
											r0, dr, dt, ntimes, wdf, dfm)
		implicit none

		! Arguments:
		real, intent(out) :: fi			! Igniting surface fire intensity (kW / sq m)
		real, intent(out) :: ti			! Igniting surface fire residence time (s)
		real, intent(out) :: u			! Mean horizontal windspeed at top of fuelbed (m/s).
		real, intent(out) :: d			! Depth of fuelbed (m)
		real, intent(out) :: tpamb		! Ambient temperature (K)
		real, intent(out) :: ak			! Area influence factor (ak / K_a parameter)
		real, intent(out) :: r0			! Minimum value of mixing parameter
		real, intent(out) :: dr			! Max - min value of mixing parameter
		real, intent(out) :: dt			! Time step for integration of burning rates (s)
		integer, intent(out):: ntimes	! Number time steps.
		real, intent(out) :: wdf		! Duff loading (kg/m^2, aka W sub d)
		real, intent(out) :: dfm		! Ratio of moisture mass to dry organic mass /
										! duff fractional moisture (aka R sub M)

		! Locals:
		real :: tamb ! Ambient temperature intermediate in celsius.
		integer :: readStat ! IO error status.

		! Constants:
		real, parameter :: fir1 = 100.0
		real, parameter :: fir2 = 1.0e+05
		real, parameter :: ti1 = 10.0
		real, parameter :: ti2 = 200.0
		real, parameter :: u1 = 0.0
		real, parameter :: u2 = 5.0
		real, parameter :: d1 = 0.1
		real, parameter :: d2 = 2.0
		real, parameter :: tam1 = -40.0
		real, parameter :: tam2 = 40.0
		real, parameter :: wdf1 = 0.1
		real, parameter :: wdf2 = 30.0
		real, parameter :: dfm1 = 0.1
		real, parameter :: dfm2 = 1.972

		fi = AskForReal2('Enter igniting surface fire intensity , kW / sq m   ', fir1, fir2)

		ti = AskForReal2('Enter igniting surface fire residence time , s   ', ti1, ti2)

		u = AskForReal2('Windspeed at top of fuelbed , m / s   ', u1, u2)

		d = AskForReal2('Depth of fuelbed , m   ', d1, d2)

		tamb = AskForReal2(' Ambient temperature , deg C   ', tam1, tam2)
		tpamb = tamb + 273.0

		do
			write(*, "(' Dimensionless area influence factor [ak parameter]   ',$)")
		
			read(*, *, iostat = readStat) ak
			if (readStat .eq. 0) then
				exit
			end if
		end do

		do
			write(*, "(' Fire environment minimum temperature parameter r0   ',$)")
		
			read(*, *, iostat = readStat) r0
			if (readStat .eq. 0) then
				exit
			end if
		end do

		do
			write(*, "(' Fire environment increment temperature parameter dr  ',$)")
	
			read(*, *, iostat = readStat) dr
			if (readStat .eq. 0) then
				exit
			end if
		end do

		do
			write(*, "(' Time step for integration of burning rates , s   ',$)")
		
			read(*, *, iostat = readStat) dt
			if (readStat .eq. 0) then
				exit
			end if
		end do

		do
			write(*, "(' Number time steps [ must exceed 1 ]   ',$)")
		
			read(*, *, iostat = readStat) ntimes
			if (readStat .eq. 0) then
				if (ntimes .le. 1) then
					cycle
				else
					exit
				end if
			end if
		end do

		wdf = AskForReal2('Duff dry weight loading, kg / sq m   ', wdf1, wdf2)

		dfm = AskForReal2('Duff moisture fraction   ', dfm1, dfm2)

	end subroutine GetFireAndSimProperties


	! Give the user the user the opportunity to review and revise data:
	!
	! History: This routine reproduces code originally in GETDAT() (block starting at label 50).
	! Testing status: Menu reproduces original behavior.
	subroutine ReviewDataMenu(parts, wdry, ash, htval, fmois, dendry, &
								sigma, cheat, condry, tpig, tchar, maxno, number, &
								fi, ti, u, d, tpamb, &
								ak, r0, dr, dt, wdf, dfm, ntimes)
		implicit none

		! Arguments:
		! 	All arguments pass through this function.
		! Fuel component properties:
		character*12, intent(inout) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(inout) :: wdry(maxno)		! Ovendry mass loading, kg/sq m
		real*4, intent(inout) :: ash(maxno)			! Mineral content, fraction dry mass
		real*4, intent(inout) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(inout) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(inout) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(inout) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(inout) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(inout) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(inout) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(inout) :: tchar(maxno)		! Char temperature, K
		integer, intent(in) :: maxno				! The maximum number of fuel classes allowed.
		integer, intent(inout) :: number			! The number of fuel classes.  May change with
													! with this call.
		! Igniting fire and environmental data:
		real*4, intent(inout) :: fi		! Current fire intensity (site avg), kW / sq m
		real*4, intent(inout) :: ti		! Igniting fire residence time (s).
		real*4, intent(inout) :: u		! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4, intent(inout) :: d		! Fuelbed depth (m)
		real*4, intent(inout) :: tpamb	! Ambient temperature (K)
		! Internal and control variables:
		real*4, intent(inout) :: ak			! Area influence factor (ak / K_a parameter)
		real*4, intent(inout) :: r0			! Minimum value of mixing parameter
		real*4, intent(inout) :: dr			! Max - min value of mixing parameter
		real*4, intent(inout) :: dt			! Time step for integration of burning rates (s)
		real*4, intent(inout) :: wdf		! Duff loading (kg/m^2, aka W sub d)
		real*4, intent(inout) :: dfm		! Ratio of moisture mass to dry organic mass /
											! duff fractional moisture (aka R sub M).
		integer, intent(inout) :: ntimes	! Number of time steps.

		! Locals:
		integer :: next ! User menu selection.
		integer :: readStat ! IO error status.

		! Constants: NA
	
		do
			! Menu:
			write(*, *) 'All data can be reviewed and revised now; Enter:'
			write(*, *) '0 to archive and/or execute present data'
			write(*, *) '1 to review fuel component properties'
			write(*, *) '2 to review igniting fire and environmental data'
			write(*, *) '3 to review internal and control variables'
			write(*, "(' 4 to terminate program now   ',$)")
		
			read(*, *, iostat = readStat) next
			if (readStat .ne. 0) then
				cycle ! If there is an error go back to the menu.
			else
			
				if (next .EQ. 0) then
					if (ArchiveMenu(parts, wdry, ash, htval, fmois, dendry, &
								sigma, cheat, condry, tpig, tchar, maxno, number, &
								fi, ti, u, d, tpamb, &
								ak, r0, dr, dt, wdf, dfm, ntimes)) then
						exit ! Return to main program.
					end if
				else if (next .EQ. 1) then
					call ReviewFuelComp(parts, wdry, ash, htval, fmois, dendry, &
										sigma, cheat, condry, tpig, tchar, maxno, number)
				else if (next .EQ. 2) then
					call ReviewFireEnvData(fi, ti, u, d, tpamb)
				else if (next .EQ. 3) then
					! Review internal and control variables:
					call ReviewIntlCntlVars(ak, r0, dr, dt, wdf, dfm, ntimes)
				else if (next .EQ. 4) then
					stop ' Terminated'
				!else ! = if ((next .LT. 1) .OR. (next .GT. 3))
				!	cycle ! Invalid value, return to menu.
				end if
			end if
		end do

	end subroutine ReviewDataMenu


	! Present a menu to allow the user to add or delete fuel types:
	!
	! History: This routine reproduces code originally in GETDAT() (block starting at label 66).
	! Note: Due to the way this code is implemented junk data may occupy invalid fuel indexes
	! after this call.  It would be better to set such data members fo an 'empty' value.
	! Deleting fuel components may result the order of elements changing.  This may not pretty,
	! but doesn't effect results because the fuels will be sorted when calculations are performed.
	! Testing status: Reproduces original behavior.
	subroutine AddDeleteComponent(parts, &
									wdry, ash, htval, fmois, dendry, &
									sigma, cheat, condry, tpig, tchar, maxno, number)
		implicit none

		! Arguments:
		character*12, intent(inout) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(inout) :: wdry(maxno)		! Ovendry mass loading, kg/sq m
		real*4, intent(inout) :: ash(maxno)			! Mineral content, fraction dry mass
		real*4, intent(inout) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(inout) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(inout) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(inout) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(inout) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(inout) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(inout) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(inout) :: tchar(maxno)		! Char temperature, K
		! Order?:
		integer, intent(in) :: maxno				! The actual number of fuel classes.
		integer, intent(inout) :: number			! The current number of fuel classes.  The
													! number of classes will be updated and returned
													! if changed by this function.

		! Locals:
		integer :: ido, ind ! User input values.
		integer :: n ! Counter
		integer :: readStat ! IO error status.

		do
			! Menu: (labeled 66 in original code)
			write(*, *) 'Enter - 1 to delete a fuel component'
			write(*, *) 'Enter 0 if no more additions or deletions'
			write(*, "(' Enter + 1 to add a fuel component   ', $)")
		
			read(*, *, iostat = readStat) ido
			if (readStat .ne. 0) then
				cycle
			end if
		
			if (ido .EQ. 0) then
				exit ! Return to the parent menu.
				! The original code effectively takes us to ReviewDataMenu().
				! Returning will take us to ReviewFuelComp(), and ReviewFuelComp() has to return
				! us to ReviewDataMenu().
			else if (ido .EQ. -1) then ! Delete a fuel component:

				do ! Menu: (labeled 69 in original code)
					! Show an abbreviated table of the fuel types:
					write(*, "(i3,3x,a12,'   loading = 'e12.3'   kg / sq m')") &
								(n, parts(n), wdry(n), n = 1, number) ! Implied do loop.
 
					write(*, "(' Index number of component to delete   ', $)")
					read(*, *, iostat = readStat) ind
					if (readStat .eq. 0) then
						exit
					end if
					! Else there is a read error, ask again.
				end do
 
				if ((ind .LE. 0) .OR. (ind .GT. number)) then ! Invalid value, start over:
					cycle ! Back to top menu.
				else if (ind .EQ. number) then
					! If the item to be removed is the last one we can just ignore that element:
					! Note: This means the data remains.
					number = number - 1
					cycle ! Back to top menu.
				else
					! The else in the orignal is implied.
				
					! Otherwise overwrite the item to be removed with the last element:
					! This changes the order but this doesn't matter because th list will be sorted
					! later.
					parts(ind) = parts(number)
					wdry(ind) = wdry(number)
					ash(ind) = ash(number)
					htval(ind) = htval(number)
					fmois(ind) = fmois(number)
					dendry(ind) = dendry(number)
					sigma(ind) = sigma(number)
					cheat(ind) = cheat(number)
					condry(ind) = condry(number)
					tpig(ind) = tpig(number)
					tchar(ind) = tchar(number)
					number = number - 1 ! JMR_NOTE: Could be combined with previous else.
					cycle ! Back to top menu.
				end if

			else if (ido .EQ. 1) then ! Add component to the end of the list:
				number = number + 1 ! Increment first to add to the end of the list.
				call GetComponentParameters(parts, wdry, ash, htval, fmois, dendry, &
											sigma, cheat, condry, tpig, tchar, number)
			
				exit ! Return to ReviewDataMenu() (see note above).
			else
				cycle ! Invalid value entered, start again:
			end if ! ido check
			! All conditions are caught.
		end do

	end subroutine AddDeleteComponent


	! Present a menu to review fuel component properties:
	!
	! History: This routine reproduces code originally in GETDAT() (block starting at label 53).
	! Testing status: Reproduces original behavior.
	subroutine ReviewFuelComp(parts, wdry, ash, htval, fmois, dendry, &
								sigma, cheat, condry, tpig, tchar, maxno, number)
		implicit none

		! Arguments:
		character*12, intent(inout) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(inout) :: wdry(maxno)		! Ovendry mass loading, kg/sq m
		real*4, intent(inout) :: ash(maxno)			! Mineral content, fraction dry mass
		real*4, intent(inout) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(inout) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(inout) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(inout) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(inout) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(inout) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(inout) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(inout) :: tchar(maxno)		! Char temperature, K
		integer, intent(in) :: maxno				! The maximum number of fuel classes allowed.
		integer, intent(inout) :: number			! The current number of fuel classes.  This may
													! change if the user opts to edit them.

		! Locals:
		integer :: readStat ! IO error status.
		integer :: ido ! User input value.
	
		do
			write(*, *) 'Enter 0 to terminate review & revision of fuel data'
			write(*, *) 'Enter 1 to select a fuel component for review'
			write(*, "(' Enter 2 to add or delete a fuel component   ',$)")

			read(*, *, iostat = readStat) ido
			if (readStat .ne. 0) then
				cycle ! If there is an error, start again.
			else
				if ((ido .lt. 0) .or. (ido .gt. 2)) then
					cycle ! Invalid value, return to menu.
				else
				
					if (ido .eq. 0) then
						exit ! Drop back to calling menu.
					else if (ido .eq. 1) then
						call ReviseFuelComponent(parts, &
													wdry, ash, htval, fmois, dendry, &
													sigma, cheat, condry, tpig, tchar, maxno, number)
					else if (ido .eq. 2) then
						call AddDeleteComponent(parts, &
													wdry, ash, htval, fmois, dendry, &
													sigma, cheat, condry, tpig, tchar, maxno, number)
						exit ! Return to ReviewDataMenu() to preserve original behavior.
					else
						cycle ! Invalid value, return to menu.
					end if
				end if
			end if
		end do

	end subroutine ReviewFuelComp


	! Allow the user to revise fuel components:
	!
	! History: This routine reproduces code originally in GETDAT() (block starting at label 55).
	! Testing status: Reproduces original behavior.
	! FuelType is shorter than FuelComponent and perhaps clearer?
	subroutine ReviseFuelComponent(parts, wdry, ash, htval, fmois, dendry, &
									sigma, cheat, condry, tpig, tchar, maxno, number)
		implicit none

		! Arguments:
		character*12, intent(inout) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(inout) :: wdry(maxno)		! Ovendry mass loading, kg/sq m
		real*4, intent(inout) :: ash(maxno)			! Mineral content, fraction dry mass
		real*4, intent(inout) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(inout) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(inout) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(inout) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(inout) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(inout) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(inout) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(inout) :: tchar(maxno)		! Char temperature, K
		integer, intent(in) :: maxno				! The maximum number of fuel classes allowed.
		integer, intent(in) :: number				! The actual number of fuel classes.

		! Locals:
		integer n ! Implied do loop.
		integer :: nn ! User selected component number.
		integer :: readStat ! IO error status.

		! Constants: NA

		do
	
			! List the fuel components:
			write(*, "(i3,3x,a12,'   loading = 'e12.3'   kg / sq m')") &
					(n, parts(n), wdry(n), n = 1, number) ! Implied do loop.

			do
				! Select the component:
				write(*, "(' Number component to be revised [ 0 = prev menu ]   ',$)")

				read(*, *, iostat = readStat) nn
				if (readStat .ne. 0) then
					cycle ! If there is an error, request number again..
				else
					if (nn .eq. 0) then
						return  ! Drop back to calling menu.
					else if ((nn .LT. 1) .or. (nn .GT. number)) then
						exit ! Invalid value, redisplay component list and menu.
					else
						call ReviseFuelParameters(parts, wdry, ash, htval, fmois, dendry, &
													sigma, cheat, condry, tpig, tchar, nn)
						exit ! On return redisplay component list and menu.
					end if
				end if
			end do
		end do ! List
	
	end subroutine ReviseFuelComponent


	! Show the user a menu of the parameters for a given fuel type (index = nFuel) and allow them
	! to edit the values.
	!
	! History: This routine reproduces code originally in GETDAT() (block starting at label 50).
	! Testing status: Reproduces original behavior.
	subroutine ReviseFuelParameters(parts, wdry, ash, htval, fmois, dendry, &
									sigma, cheat, condry, tpig, tchar, nFuel)
		implicit none

		! Arguments:
		character*12, intent(inout) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(inout) :: wdry(maxno)		! Ovendry mass loading, kg/sq m
		real*4, intent(inout) :: ash(maxno)			! Mineral content, fraction dry mass
		real*4, intent(inout) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(inout) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(inout) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(inout) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(inout) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(inout) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(inout) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(inout) :: tchar(maxno)		! Char temperature, K
		integer, intent(in) :: nFuel				! The fuel component to modify.

		! Locals:
		real :: idn, value ! User input values.
		integer :: readStat ! IO error status.

		! Constants: NA
	
		do
			! Menu of parameters:
			write(*, "(' index parameter data entered')")
			write(*, "('    1'10x'name'4x,a12)") parts(nFuel)
			write(*, "('    2'7x'loading  ',e12.3)") wdry(nFuel)
			write(*, "('    3'7x'ash frac ',e12.3)") ash(nFuel)
			write(*, "('    4'7x'heat comb',e12.3)") htval(nFuel)
			write(*, "('    5'7x'mois frac',e12.3)") fmois(nFuel)
			write(*, "('    6'7x'density  ',e12.3)") dendry(nFuel)
			write(*, "('    7'7x'surf/vol ',e12.3)") sigma(nFuel)
			write(*, "('    8'7x'heat capy',e12.3)") cheat(nFuel)
			write(*, "('    9'7x'conduct y',e12.3)") condry(nFuel)
			write(*, "('   10'7x'ig temp C',e12.3)") (tpig(nFuel) - 273.0)
			write(*, "('   11'7x'char temp',e12.3)") (tchar(nFuel) - 273.0)
			write(*, "('   enter number param to change [ 0 = prev menu ]   ',$)")
	
			read(*, *, iostat = readStat) idn
		
			if (readStat .ne. 0) then
				cycle ! If there is an error, start again.
			else
				if (idn .eq. 0) then
					exit ! (= return) Return to parent menu (55)
				else if (idn .GT. 11) then
					cycle ! Invalid value, return to menu. (59)
				else if (idn .EQ. 1) then
					write(*, *) 'Enter new name - 12 char or fewer'
					read(*, "(a12)") parts(nFuel)
					!cycle ! Allow more properties to be edited.
				else ! Items 2-11:
				
					!write(*, "format(' New value = ?   '$)")
					write(*, "(' New value = ?   '$)")
				
					read(*, *, iostat = readStat) value
					if (readStat .ne. 0) then
						cycle ! If there is an error, start again.
					else
						! Store the revised value:
						if (idn .EQ. 2) then
							wdry(nFuel) = value
						else if (idn .EQ. 3) then
							ash(nFuel) = value
						else if (idn .EQ. 4) then
							htval(nFuel) = value
						else if (idn .EQ. 5) then
							fmois(nFuel) = value
						else if (idn .EQ. 6) then
							dendry(nFuel) = value
						else if (idn .EQ. 7) then
							sigma(nFuel) = value
						else if (idn .EQ. 8) then
							cheat(nFuel) = value
						else if (idn .EQ. 9) then
							condry(nFuel) = value
						else if (idn .EQ. 10) then
							tpig(nFuel) = value + 273.0
						else if (idn .EQ. 11) then
							tchar(nFuel) = value + 273.0
						end if
					end if
				end if
			end if
		end do

	end subroutine ReviseFuelParameters


	! Present a menu to review igniting fire and environmental data:
	!
	! History: This routine reproduces code originally in GETDAT() (block starting at label 71).
	! Testing status: Incomplete.
	subroutine ReviewFireEnvData(fi, ti, u, d, tpamb)
		implicit none

		! Arguments:
		real*4, intent(inout) :: fi		! Current fire intensity (site avg), kW / sq m
		real*4, intent(inout) :: ti		! Igniting fire residence time (s).
		real*4, intent(inout) :: u		! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4, intent(inout) :: d		! Fuelbed depth (m)
		real*4, intent(inout) :: tpamb	! Ambient temperature (K)

		! Locals:
		integer :: readStat	! IO error status.
		integer :: ind		! User menu input selection.
		real*4 :: value	! User input value.

		! Constants: NA

		do
			! Parameter selection menu:
			write(*, "('  1 - Igniting fire intensity, kW/ sq m      'e12.3)") fi
			write(*, "('  2 - Igniting fire residence time, s        'e12.3)") ti
			write(*, "('  3 - windspeed at top of fuel bed, m / s    'e12.3)") u
			write(*, "('  4 - fuelbed depth, m                       'e12.3)") d
			write(*, "('  5 - ambient temperature, deg C             'e12.3)") (tpamb - 273.0)

			do
				write(*, "(' Index of parameter to change [ 0 = prev menu ]   ',$)")

				! Read the user's selection:
				read(*, *, iostat = readStat) ind 
				if (readStat .eq. 0) then
					exit ! Success, proceed to next step.
				end if
				! If there was an error prompt again.
			end do

			! Respond:
			if (ind .EQ. 0) then
				return ! Return to calling menu.
			else if ((ind .LT. 0) .OR. (ind .GT. 5)) then
				cycle ! For invalid menu selections redisplay the menu.
			else
				! Get and save the new value:
				write(*, "(' New value = ?   '$)")
				
				read(*, *, iostat = readStat) value
				if (readStat .ne. 0) then
					cycle ! If there is an error, start again.
				else
					! Store the revised value:
					if (ind .EQ. 1) then
						fi = value
					else if (ind .EQ. 2) then
						ti = value
					else if (ind .EQ. 3) then
						u = value
					else if (ind .EQ. 4) then
						d = value
					else if (ind .EQ. 5) then
						tpamb = value + 273.0
					end if
				end if
				! Now return to the menu for further changees.
			end if
		end do

	end subroutine ReviewFireEnvData


	! Present a menu to review internal and control variables:
	!
	! History: This routine reproduces code originally in GETDAT() (block starting at label 79).
	! Testing status: Incomplete.
	subroutine ReviewIntlCntlVars(ak, r0, dr, dt, wdf, dfm, ntimes)
		implicit none

		! Arguments:
		real*4, intent(inout) :: ak			! Area influence factor (ak / K_a parameter)
		real*4, intent(inout) :: r0			! Minimum value of mixing parameter
		real*4, intent(inout) :: dr			! Max - min value of mixing parameter
		real*4, intent(inout) :: dt			! Time step for integration of burning rates (s)
		real*4, intent(inout) :: wdf		! Duff loading (kg/m^2, aka W sub d)
		real*4, intent(inout) :: dfm		! Ratio of moisture mass to dry organic mass /
											! duff fractional moisture (aka R sub M)
		integer, intent(inout) :: ntimes	! Number of time steps.

		! Locals:
		integer :: readStat	! IO error status.
		integer :: ind		! User menu input selection.
		real*4 :: value	! User input value.

		! Constants: NA

		do
			! Parameter selection menu:
			write(*, "('  1 - Area influence factor [ ak parameter ]       'e12.3)") ak
			write(*, "('  2 - Fire environment temperature parameter r0    'e12.3)") r0
			write(*, "('  3 - Fire environment temperature parameter dr    'e12.3)") dr
			write(*, "('  4 - Time step to integrate burning rates , s     'e12.3)") dt
			write(*, "('  5 - Duff dry weight loading , kg / sq m          'e12.3)") wdf
			write(*, "('  6 - Duff moisture content , fraction dry weight  'e12.3)") dfm
			write(*, "('  7 - Number time steps                            'i4)") ntimes
			write(*, "('  Index of parameter to change [ 0 = prev menu ]   ',$)")

			! Read the user's selection:
			read(*, *, iostat = readStat) ind ! value
			if (readStat .ne. 0) then
				cycle
			end if
		
			! Respond:
			if (ind .EQ. 0) then
				exit ! Return to calling menu.
			else if ((ind .LT. 0) .OR. (ind .GT. 7)) then
				cycle ! For invalid menu selections redisplay the menu.
			else
				! Get and save the new value:
				write(*, "(' New value = ?   '$)")
			
				if (ind .LT. 7) then
					read(*, *, iostat = readStat) value
					if (readStat .ne. 0) then
						cycle ! If there is an error, start again.
					else
						! Store the revised value:
						if (ind .EQ. 1) then
							ak = value
						else if (ind .EQ. 2) then
							r0 = value
						else if (ind .EQ. 3) then
							dr = value
						else if (ind .EQ. 4) then
							dt = value
						else if (ind .EQ. 5) then
							wdf = value
						else if (ind .EQ. 6) then
							dfm = value
						end if
					end if
				else ! ind == 7
					read(*, *, iostat = readStat) ntimes
					if (readStat .ne. 0) then
						cycle ! If there is an error, start again.
					end if
				end if
			end if
		end do

	end subroutine ReviewIntlCntlVars


	! Retrive data from input files (stored previously with ArchiveSettings():
	!
	! History: This routine is based on code originally in GETDAT() (block starting at label 1000).
	! Testing status: Reproduces original behavior.
	subroutine GetDataFromFiles(parts, wdry, ash, htval, &
									fmois, dendry, sigma, cheat, condry, tpig, &
									tchar, number, maxno, &
									fi, ti, u, d, tpamb, &
									ak, r0, dr, dt, ntimes, wdf, dfm)
		implicit none

		! Arguments:
		! Match order to other functions!
		character*12, intent(out) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(out) :: wdry(maxno) 		! Ovendry mass loading, kg/sq m
		real*4, intent(out) :: ash(maxno)		! Mineral content, fraction dry mass
		real*4, intent(out) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(out) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(out) :: dendry(maxno)	! Ovendry mass density, kg / cu m
		real*4, intent(out) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(out) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(out) :: condry(maxno)	! Thermal conductivity, W / m K, ovendry
		real*4, intent(out) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(out) :: tchar(maxno)		! Char temperature, K
		integer, intent(inout) :: number		! The actual number of fuel classes
		integer, intent(in) :: maxno			! The maximum number of fuel classes allowed.
		real*4, intent(out) :: fi				! Current fire intensity (site avg), kW / sq m
		real*4, intent(out) :: ti				! Igniting fire residence time (s)
		real*4, intent(out) :: u				! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4, intent(out) :: d				! Fuelbed depth (m)
		real*4, intent(out) :: tpamb			! Ambient temperature (K)
	
		real*4, intent(out) :: ak				! Area influence factor (ak / K_a parameter)
		real*4, intent(out) :: r0				! Minimum value of mixing parameter
		real*4, intent(out) :: dr				! Max - min value of mixing parameter
		real*4, intent(out) :: dt				! Time step for integration of burning rates (s)
		integer, intent(out) :: ntimes			! Number of time steps.
		real*4, intent(out) :: wdf				! Duff loading (kg/m^2, aka W sub d)
		real*4, intent(out) :: dfm				! Ratio of moisture mass to dry organic mass /
												! duff fractional moisture (aka R sub M)

		! Locals:
		integer :: readStat, openStat			! IO error statuses.
		character*12 :: infile					! Stores the name of input data files.
		integer :: n							! Counter.

		! Constants:
		character(len = *), parameter :: format08 = "(a12)"
		character(len = *), parameter :: format2004 = "(2x,a12)"
		character(len = *), parameter :: format2005 = "(2x,5e15.3/2x,5e15.3)"
		character(len = *), parameter :: format2008 = "(2x,5e15.3)"
		character(len = *), parameter :: format2009 = "(2x,4e15.3, i6)"
		character(len = *), parameter :: format2010 = "(2x,2e15.3)"
		integer, parameter :: fun = 77 ! File unit number.  Reused for both files.

		! Prompt for the first file:
		do
			write(*, "(' Enter name of file holding fuel component data   ',$)")

			read(*, format08, iostat = readStat) infile
			if (readStat .eq. 0) then
				exit
			end if
		end do

		! Open the file:
		do
			open(fun, file = infile, status = 'OLD', form= 'FORMATTED', iostat = openStat)
			if (openStat .eq. 0) then
				exit
			else
				do
					write(*, "(' Error opening input file  ', a12'  reenter file name')") infile
					read(*, format08, iostat = readStat) infile
					if (readStat .eq. 0) then
						exit ! Exit to out loop and try to open the new file.
					end if
				end do
			end if
		end do

		! Read in the data:
		read(fun, "(i6)") number
		do n = 1, number
			read(fun, format2004) parts(n)
			read(fun, format2005) wdry(n), ash(n), htval(n), fmois(n), &
									dendry(n), sigma(n), cheat(n), &
									condry(n), tpig(n), tchar(n)
		end do
		close(fun)


		! Prompt for the second file:
		! This is so ugly!
		write(*, "(' Enter name of file holding igniting fire, environmental'/"// &
				"', and program control parameters   ',$)")
	
		read(*, format08) infile
		! JMR_Note: No error checking for this file!

		! Open the file:
		do
			open(fun, file = infile, status = 'OLD', form= 'FORMATTED', iostat = openStat)
			if (openStat .eq. 0) then
				exit ! Opened successfully.
			else
				do
					write(*, "(' Error opening input file  ', a12'  reenter file name')") infile
					read(*, format08, iostat = readStat) infile
					if (readStat .eq. 0) then
						exit ! Exit to out loop and try to open the new file.
					end if
				end do
			end if
		end do
	
		! Read in the data:
		read(fun, format2008) fi, ti, u, d, tpamb
		read(fun, format2009) ak, r0, dr, dt, ntimes
		read(fun, format2010, iostat = readStat) wdf, dfm
		! JMR_NOTE: I'm not sure why these two values may be missing.
		!if (readStat .eq. iostat_end) then
		! JMR_NOTE: The implementation of iostat_end may be compiler specific.  The following check
		! should be valid for GFortran but may not be for others.
		if (is_iostat_end(readStat)) then
			wdf = 0.0
			dfm = 0.0
		end if

		close(fun)

	end subroutine GetDataFromFiles


	! Prompt the user about storing the current settings to file:
	!
	! History: This routine reproduces code originally in GETDAT() (block starting at label 2000).
	! Testing status: Reproduces original behavior.
	function ArchiveMenu(parts, wdry, ash, htval, fmois, dendry, &
								sigma, cheat, condry, tpig, tchar, maxno, number, &
								fi, ti, u, d, tpamb, &
								ak, r0, dr, dt, wdf, dfm, ntimes) result(exitMenus)

		implicit none
	
		! Arguments: All argument pass through to ArchiveSettings().
		character*12, intent(in) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(in) :: wdry(maxno) 		! Ovendry mass loading, kg/sq m
		real*4, intent(in) :: ash(maxno)		! Mineral content, fraction dry mass
		real*4, intent(in) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(in) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(in) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(in) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(in) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(in) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(in) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(in) :: tchar(maxno)		! Char temperature, K
		integer, intent(in) :: maxno			! The maximum number of fuel classes allowed.
		integer, intent(in) :: number			! The actual number of fuel classes.
		real*4, intent(in) :: fi				! Current fire intensity (site avg), kW / sq m
		real*4, intent(in) :: ti				! Igniting fire residence time (s)
		real*4, intent(in) :: u					! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4, intent(in) :: d					! Fuelbed depth (m)
		real*4, intent(in) :: tpamb				! Ambient temperature (K)
		real*4, intent(in) :: ak				! Area influence factor (ak / K_a parameter)
		real*4, intent(in) :: r0				! Minimum value of mixing parameter
		real*4, intent(in) :: dr				! Max - min value of mixing parameter
		real*4, intent(in) :: dt				! Time step for integration of burning rates (s)
		real*4, intent(in) :: wdf				! Duff loading (kg/m^2, aka W sub d)
		real*4, intent(in) :: dfm				! Ratio of moisture mass to dry organic mass /
												! duff fractional moisture (aka R sub M)
		integer, intent(in) :: ntimes			! Number of time steps.

		! Locals:
		logical :: exitMenus ! Should the parent menus exit to the base of the menus?
		integer :: readStat ! IO error status.
		integer :: ido ! User menu input selection.

		! Constants:= NA

		exitMenus = .false.

		do
			write(*, *) 'Enter 2 for previous menu'
			write(*, *) 'Enter 1 to archive current data set'
			write(*, "(' Enter 0 to execute without archiving   '$)")
		
			read(*, *, iostat = readStat) ido
			if (readStat .eq. 0) then
				if (ido .eq. 0) then
					exitMenus = .true.
					return
				else if (ido .eq. 1) then
					call ArchiveSettings(parts, wdry, ash, htval, fmois, dendry, &
												sigma, cheat, condry, tpig, tchar, maxno, number, &
												fi, ti, u, d, tpamb, &
												ak, r0, dr, dt, wdf, dfm, ntimes)
				else if (ido .eq. 2) then
					return ! Drop back to calling menu (ReviewDataMenu()).
				!else
				!	cycle ! Invalid menu value.
				end if
			end if
		end do

	end function ArchiveMenu


	! Save the current run conditions (settings and parameters) to file:
	!
	! This code was pulled out of ArchiveMenu() to increase readability.
	! Testing status: Produces identical output files to original code.
	subroutine ArchiveSettings(parts, wdry, ash, htval, fmois, dendry, &
								sigma, cheat, condry, tpig, tchar, maxno, number, &
								fi, ti, u, d, tpamb, &
								ak, r0, dr, dt, wdf, dfm, ntimes) ! Name SaveSettings?

		implicit none
	
		! Arguments:
		! Reorder to match parent routines?
		character*12, intent(in) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(in) :: wdry(maxno) 		! Ovendry mass loading, kg/sq m
		real*4, intent(in) :: ash(maxno)		! Mineral content, fraction dry mass
		real*4, intent(in) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(in) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(in) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(in) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(in) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(in) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(in) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(in) :: tchar(maxno)		! Char temperature, K
		integer, intent(in) :: maxno			! The maximum number of fuel classes allowed.
		integer, intent(in) :: number			! The actual number of fuel classes.
		real*4, intent(in) :: fi				! Current fire intensity (site avg), kW / sq m
		real*4, intent(in) :: ti				! Igniting fire residence time (s)
		real*4, intent(in) :: u					! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4, intent(in) :: d					! Fuelbed depth (m)
		real*4, intent(in) :: tpamb				! Ambient temperature (K)
		real*4, intent(in) :: ak				! Area influence factor (ak / K_a parameter)
		real*4, intent(in) :: r0				! Minimum value of mixing parameter
		real*4, intent(in) :: dr				! Max - min value of mixing parameter
		real*4, intent(in) :: dt				! Time step for integration of burning rates (s)
		real*4, intent(in) :: wdf				! Duff loading (kg/m^2, aka W sub d)
		real*4, intent(in) :: dfm				! Ratio of moisture mass to dry organic mass /
												! duff fractional moisture (aka R sub M)
		integer, intent(in) :: ntimes			! Number of time steps.

		! Locals:
		integer :: readStat, openstat	! IO error statuses.
		character*12 :: outfil			! Stores the name of input data files.
		integer :: n					! Counter.
		integer :: ido					! User menu input selection.
	
		! Constants:
		character(len = *), parameter :: format08 = "(a12)"
		character(len = *), parameter :: format2004 = "(2x,a12)"
		character(len = *), parameter :: format2005 = "(2x,5e15.3/2x,5e15.3)"
		character(len = *), parameter :: format2008 = "(2x,5e15.3)"
		character(len = *), parameter :: format2009 = "(2x,4e15.3, i6)"
		character(len = *), parameter :: format2010 = "(2x,2e15.3)"
		integer, parameter :: fun = 66 ! File unit number.  Reused for both files.

		! Archiving:
		do
			write(*, *) 'Enter file name [ 12 char or fewer ] for storage of the'
			write(*, "(' fuel component data currently in use   ',$)")
		
			read(*, format08, iostat = readStat) outfil
			if (readStat .eq. 0) then
				open(fun, file = outfil, status = 'UNKNOWN', form = 'FORMATTED', iostat = openStat)
				if (openStat .ne. 0) then
					write(*, "(' Error opening output file  'a12'  enter another name')") outfil
				else
					exit ! File opened, proceed to writing.
				end if
			end if
		end do

		! Write to fuel component file:
		write(fun, "(i6,' Fuel component properties -- 3 records per entry')") number

		do n = 1, number
			write(fun, format2004) parts(n)
			write(fun, format2005) wdry(n), ash(n), htval(n), fmois(n), &
									dendry(n), sigma(n), cheat(n), &
									condry(n), tpig(n), tchar(n)
		end do

		close(fun)

		! Open the IFEC file:
	
		do
			write(*, *) 'Enter file name [ 12 char or fewer ] for storage of the'
			write(*, "(' igniting fire, environmental, and control data used  ',$)")

			read(*, format08, iostat = readStat) outfil
			if (readStat .eq. 0) then
				open(fun, file = outfil, status = 'UNKNOWN', form = 'FORMATTED', iostat = openStat)
				if (openStat .ne. 0) then
					do
						write(*, "(' Error opening output file  'a12'  enter another name')") outfil
				
						read(*, format08, iostat = readStat) outfil
						if (readStat .eq. 0) then
							exit ! Drop back and try to open the new file name.
						end if
					end do
				else
					exit ! File opened, proceed to writing.
				end if
			end if
		end do
	
		! Write data:
		write(fun, format2008) fi, ti, u, d, tpamb
		write(fun, format2009) ak, r0, dr, dt, ntimes
		write(fun, format2010) wdf, dfm
	
		close(fun)

	end subroutine ArchiveSettings


	!
	subroutine SaveFuelsToFile(parts, wdry, ash, htval, fmois, dendry, &
								sigma, cheat, condry, tpig, tchar, number)!maxno, number)!, &
								!fi, ti, u, d, tpamb, &
								!ak, r0, dr, dt, wdf, dfm, ntimes) ! Name SaveSettings?

		implicit none
	
		! Arguments:
		! Reorder to match parent routines?
		character*12, intent(in) :: parts(maxno)	! Fuel component names / labels
		real*4, intent(in) :: wdry(maxno) 		! Ovendry mass loading, kg/sq m
		real*4, intent(in) :: ash(maxno)		! Mineral content, fraction dry mass
		real*4, intent(in) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(in) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(in) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(in) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(in) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(in) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(in) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(in) :: tchar(maxno)		! Char temperature, K
		!integer, intent(in) :: maxno			! The maximum number of fuel classes allowed.
 		integer, intent(in) :: number			! The actual number of fuel classes.
! 		real*4, intent(in) :: fi				! Current fire intensity (site avg), kW / sq m
! 		real*4, intent(in) :: ti				! Igniting fire residence time (s)
! 		real*4, intent(in) :: u					! Mean horizontal windspeed at top of fuelbed (m/s).
! 		real*4, intent(in) :: d					! Fuelbed depth (m)
! 		real*4, intent(in) :: tpamb				! Ambient temperature (K)
! 		real*4, intent(in) :: ak				! Area influence factor (ak / K_a parameter)
! 		real*4, intent(in) :: r0				! Minimum value of mixing parameter
! 		real*4, intent(in) :: dr				! Max - min value of mixing parameter
! 		real*4, intent(in) :: dt				! Time step for integration of burning rates (s)
! 		real*4, intent(in) :: wdf				! Duff loading (kg/m^2, aka W sub d)
! 		real*4, intent(in) :: dfm				! Ratio of moisture mass to dry organic mass /
! 												! duff fractional moisture (aka R sub M)
! 		integer, intent(in) :: ntimes			! Number of time steps.

		! Locals:
		integer :: openstat	! IO error statuses.
		character*12 :: outfil			! Stores the name of input data files.
		integer :: n					! Counter.
		!integer :: ido					! User menu input selection.
	
		! Constants:
		!character(len = *), parameter :: format08 = "(a12)"
		character(len = *), parameter :: format2004 = "(2x,a12)"
		character(len = *), parameter :: format2005 = "(2x,5e15.3/2x,5e15.3)"
		! character(len = *), parameter :: format2008 = "(2x,5e15.3)"
! 		character(len = *), parameter :: format2009 = "(2x,4e15.3, i6)"
! 		character(len = *), parameter :: format2010 = "(2x,2e15.3)"
		integer, parameter :: fun = 66 ! File unit number.  Reused for both files.
		
		outfil = "FC_Data.txt"

		! ...
		open(fun, file = outfil, status = 'UNKNOWN', form = 'FORMATTED', iostat = openStat)
		if (openStat .ne. 0) then
			write(*, "(' Error opening output file  'a12'  not saving.')") outfil
			return
		end if

		! Write to fuel component file:
		write(fun, "(i6,' Fuel component properties -- 3 records per entry')") number

		do n = 1, number
			write(fun, format2004) parts(n)
			write(fun, format2005) wdry(n), ash(n), htval(n), fmois(n), &
									dendry(n), sigma(n), cheat(n), &
									condry(n), tpig(n), tchar(n)
		end do

		close(fun)

	end subroutine SaveFuelsToFile


	!...
	subroutine SaveFuelsToFileR(wdry, ash, htval, fmois, dendry, &
								sigma, cheat, condry, tpig, tchar, number) bind(C, name = "savefuelstofiler")
		implicit none

		! Arguments:
		double precision, intent(in) :: wdry(maxno)		! Ovendry mass loading, kg/sq m
  		double precision, intent(in) :: ash(maxno)		! Mineral content, fraction dry mass
		double precision, intent(in) :: htval(maxno)	! Low heat of combustion, J / kg
		double precision, intent(in) :: fmois(maxno)	! Moisture fraction of component
		double precision, intent(in) :: dendry(maxno)	! Ovendry mass density, kg / cu m
		double precision, intent(in) :: sigma(maxno)	! Surface to volume ratio, 1 / m
		double precision, intent(in) :: cheat(maxno)	! Specific heat capacity, (J / K) / kg dry mass
		double precision, intent(in) :: condry(maxno)	! Thermal conductivity, W / m K, ovendry
		double precision, intent(in) :: tpig(maxno)		! Ignition temperature, K
		double precision, intent(in) :: tchar(maxno)	! Char temperature, K
		integer, intent(in) :: number	! The number of fuel classes.
		
		character*12 :: parts(maxno) = "Test"
		
		call SaveFuelsToFile(parts, &
								real(wdry), real(ash), real(htval), real(fmois), real(dendry), &
								real(sigma), real(cheat), real(condry), real(tpig), real(tchar), number)

	end subroutine SaveFuelsToFileR


	! ...
	! Argument order differs from ArchiveSettings() and follows Simulate().
! 	subroutine ArchiveSettingsR
! 		(fi, ti, u, d, tpamb, ak, r0, dr, dt, wdf, dfm, ntimes, number, &
! 							wdry, ash, htval, fmois, dendry, sigma, cheat, condry, tpig, tchar, &
! 							xmat, tign, tout, wo, diam) bind(C, name = "archivesettingsR")
! 		implicit none
! 
! 		! Arguments:
! 		double precision, intent(in) :: fi		! Current fire intensity (site avg), kW / sq m
! 		double precision, intent(in) :: ti			! Igniting fire residence time (s).
! 		double precision, intent(in) :: u			! Mean horizontal windspeed at top of fuelbed (m/s).
! 		double precision, intent(in) :: d			! Fuelbed depth (m)
! 		double precision, intent(in) :: tpamb		! Ambient temperature (K)
! 		double precision, intent(in) :: ak			! Area influence factor (ak / K_a parameter)
! 		double precision, intent(in) :: r0			! Minimum value of mixing parameter
! 		double precision, intent(in) :: dr			! Max - min value of mixing parameter
! 		double precision, intent(in) :: dt			! Time step for integration of burning rates (s)
! 		double precision, intent(in) :: wdf			! Duff loading (kg/m^2, aka W sub d)
! 		double precision, intent(in) :: dfm			! Ratio of moisture mass to dry organic mass /
! 													! duff fractional moisture (aka R sub M).
! 		integer, intent(in) :: ntimes						! Number of time steps.
! 		integer, intent(in) :: number	! The number of fuel classes. ! Could try to remove?????
! 		
! 		double precision, intent(in) :: wdry(maxno)		! Ovendry mass loading, kg/sq m
! 		double precision, intent(in) :: ash(maxno)		! Mineral content, fraction dry mass
! 		double precision, intent(in) :: htval(maxno)		! Low heat of combustion, J / kg
! 		double precision, intent(in) :: fmois(maxno)		! Moisture fraction of component
! 		double precision, intent(in) :: dendry(maxno)	! Ovendry mass density, kg / cu m
! 		double precision, intent(in) :: sigma(maxno)		! Surface to volume ratio, 1 / m
! 		double precision, intent(in) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
! 		double precision, intent(in) :: condry(maxno)	! Thermal conductivity, W / m K, ovendry
! 		double precision, intent(in) :: tpig(maxno)		! Ignition temperature, K
! 		double precision, intent(in) :: tchar(maxno)		! Char temperature, K
! 		
! 		call ArchiveSettings(real(), real(), real(), real(), real(), real(), real())
! 
! 	end subroutine ArchiveSettingsR


	! Check if the value passed is in the supplied value range and offer the user a chance to
	! revise it if it is not.  True is returned if the user selected to reenter the value.
	!
	! History: Modernized original Burnup subroutine.
	subroutine RETRY(value, valo, vahi, test)
		implicit none

		! Arguments:
		real, intent(inout) :: value ! Value to check and potentially revise.
		real, intent(in) :: valo ! Low value limit.
		real, intent(in) :: vahi ! High value limit.
		logical, intent(out) :: test ! Is the value (potentially revised) in range?

		! Locals:
		integer :: in ! Menu return value.
		integer :: readStat ! IO error status.

		! The value is in the expected range.  Nothing more to do.
		test = (value .GT. valo) .AND. (value .LT. vahi)
		if (test) then
			test = .FALSE.
			return
		end if

		! The value is lower than expected:
		test = (value .LE. valo)
		if (test) then
			do !while (in .NE. 1)
				write(*, *) ' Value entered seems small'
				write(*, *) ' Enter 0 to proceed with small value'
				write(*, *) ' Enter 1 to provide a new value'

				read(*, *, iostat = readStat) in
				if (readStat .eq. 0) then 
					if (in .EQ. 0) then
						test = .FALSE.
						return
					else if (in .EQ. 1) then
						return
					! else if (in .NE. 1): Repeat for invalid menu selections.
					end if
				end if
			end do
		end if

		! The value is higher than expected:
		test = (value .ge. vahi)
		if (test) then
			do
				write(*,*) ' Value entered seems large'
				write(*,*) ' Enter 0 to proceed with large value'
				write(*,*) ' Enter 1 to provide a new value'

				read(*, *, iostat = readStat) in
				if (readStat .eq. 0) then 
					if (in .EQ. 0) then
						test = .FALSE.
						return
					else if (in .EQ. 1) then
						return
					! else if (in .NE. 1): Repeat for invalid menu selections.
					end if
				end if
			end do
		end if

	end subroutine RETRY


! -- Pagebreak --
! Pg. 88:  ARRAYS() starts here.


	!c Orders the fuel description arrays according to the paradigm described in
	!c subroutine SORTER and computes the interaction matrix xmat from the array
	!c elam and the list alone returned from subroutine OVLAPS.
	!
	! History: Modernized original Burnup subroutine.
	subroutine ARRAYS(maxno, number, wdry, ash, dendry, fmois, &
						sigma, htval, cheat, condry, tpig, tchar, &
						diam, key, work, ak, elam, alone, xmat, &
						wo, maxkl, parts, list, area)
		implicit none

		! Arguments:
		integer, intent(in) :: maxno				! The maximum number of fuel classes allowed.
		integer, intent(in) :: number				! The actual number of fuel classes.
		real*4, intent(inout) :: wdry(maxno)		! Ovendry mass loading, kg/ sq m				! JMR: Transcription!!!!
		real*4, intent(inout) :: ash(maxno)			! Mineral content, fraction dry mass
		real*4, intent(inout) :: dendry(maxno)		! Ovendry mass density, kg / cu m
		real*4, intent(inout) :: fmois(maxno)		! Moisture content, fraction dry mass
		real*4, intent(inout) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(inout) :: htval(maxno)		! Low heat of combustion, J / kg
		real*4, intent(inout) :: cheat(maxno)		! Specific heat capacity, (J / K) / kg dry mass
		real*4, intent(inout) :: condry(maxno)		! Thermal conductivity, W / m K, ovendry
		real*4, intent(inout) :: tpig(maxno)		! Ignition temperature, K
		real*4, intent(inout) :: tchar(maxno)		! Char temperature, K
		real*4, intent(out) :: diam(maxkl)			! initial diameter, m [by interaction pairs]
		integer, intent(out) :: key(maxno)			! Ordered index list
		real*4, intent(out) :: work(maxno)			! Workspace array
		real*4, intent(in) :: ak					! Area influence factor [ak parameter]
		real*4, intent(out) :: elam(maxno, maxno)	! Interaction matrix from OVLAPS
		real*4, intent(out) :: alone(maxno)			! Noninteraction fraction list from OVLAPS
		real*4, intent(out) :: xmat(maxkl)			! Consolidated interaction matrix
		real*4, intent(out) :: wo(maxkl)			! Initial dry loading by interaction pairs
		integer, intent(in) :: maxkl				! Max triangular matrix size.
		character*12, intent(inout) :: parts(maxno)	! Fuel component names / labels
		character*12, intent(out) :: list(maxno)	! Intermediary for reordering parts name array?????
													! This is passed in but is not initialized prior
													! to that.  It doesn't appear that it is used
													! after it is returned.  It appears to only be
													! used internal to this routine.
		real*4, intent(inout) :: area(maxno)		! Fraction of site area expected to be covered at
													! least once by initial planform area of ea size

		! Locals:
		integer :: j, k, kl, kj ! Counters
		real*4 :: diak, wtk

		! Testing was done to confrim that explicit initialization of locals was not needed here.

		call SORTER(maxno, number, sigma, fmois, dendry, key)

		do j = 1, number
			k = key(j)
			list(j) = parts(k)
		end do
		do j = 1, number
			parts(j) = list (j)
		end do
		do j = 1, number


! -- Pagebreak --
! Pg. 89:


			k = key(j)
			work(j) = wdry(k)
		end do
		do j = 1, number
			wdry(j) = work(j)
		end do

		do j = 1, number
			k = key (j)						! JMR: Whitespace!!!!!
			work(j) = ash(k)
		end do
		do j = 1, number
			ash(j) = work(j)
		end do

		do j = 1, number
			k = key (j)						! JMR: Whitespace!!!!!
			work(j) = htval(k)
		end do
		do j = 1, number
			htval(j) = work(j)
		end do

		do j = 1, number
			k = key(j)
			work(j) = cheat(k)
		end do
		do j = 1, number
			cheat(j) = work(j)
		end do

		do j = 1, number
			k = key (j)						! JMR: Whitespace!!!!!
			work(j) = condry(k)
		end do
		do j = 1, number
			condry(j) = work(j)
		end do

		do j = 1, number
			k = key(j)
			work(j) = tpig(k)
		end do
		do j = 1, number
			tpig(j) =  work (j)
		end do

		do j = 1, number
			k = key (j)						! JMR: Whitespace!!!!!
			work(j) = tchar(k)
		end do
		do j = 1, number
			tchar(j) = work(j)
		end do

! -- Pagebreak --
! Pg. 90:


		call OVLAPS(wdry, sigma, dendry, ak, number, maxno, maxkl, &
					fmois, &
					xmat, elam, alone, area)

		do k = 1, number
			diak = 4.0 / sigma(k)
			wtk = wdry(k)
			
			! Populate the alone/no companion indexes of the arrays:
			kl = Loc(k, 0)
			diam(kl) = diak
			xmat(kl) = alone(k)
			wo(kl) = wtk * xmat(kl)
			
			! Populate the interacting indexes of the arrays:
			do j = 1, k
				kj = Loc(k, j)
				diam(kj) = diak
				xmat(kj) = elam (k, j)					! JMR: Transcription!!!!!!
				wo(kj) = wtk * xmat(kj)
			end do
		end do

		! JMR: Temporary reporting:
		! print *, "ARRAYS()"
! 		print *, "diam:"
! 		print *, diam
! 		print *, "xmat:"
! 		print *, xmat
! 		print *, "wo:"
! 		print *, wo

	end subroutine ARRAYS


! -- Pagebreak --
! Pg. 91:


	!c Sorts fuel element list in order of increasing size (decreasing sigma)
	!c For elements with same size order determined on increasing moisture
	!c content (fmois). If items have same size and moisture content, order
	!c on the basis of increasing mass density (dryden). "number" elements are
	!c included in the list, which has a maximum length of "maxno". The integer
	!c list: key(j), j = 1, number holds the indices in order, so other
	!c fuel parameters can be ordered and associated as necessary.
	!
	! History: Modernized original Burnup subroutine.
	subroutine SORTER(maxno, number, sigma, fmois, dryden, key)
		implicit none

		! Arguments:
		integer, intent(in) :: maxno				! The maximum number of fuel classes allowed.
		integer, intent(in) :: number				! The actual number of fuel classes.
		real*4, intent(inout) :: sigma(maxno)		! Surface to volume ratio, 1 / m
		real*4, intent(inout) :: fmois(maxno)		! Moisture content, fraction dry mass
		real*4, intent(inout) :: dryden(maxno)		! Ovendry mass density, kg / cu m
		integer, intent(inout) :: key(maxno)		! Ordered index list

		! Locals:
		integer :: j, i ! Counters
		real :: s, fm, de, keep, usi
		logical :: diam, mois, dens, tied
		logical :: newIndexFound ! Note: Not present in original code.

		newIndexFound = .false.

		do j = 1, maxno
			key (j) = j
		end do

		!c Replacement sort: order on increasing size, moisture, density
		do j = 2, number
			! Store the values for this fuel index:
			s = 1.0 / sigma (j)					! JMR: Whitespace!!!!!
			fm = fmois(j)
			de = dryden(j)
			keep = key(j)
	
			! Compare this index (j) with every index before it:
			do i = (j - 1), 1, -1
				usi = 1.0 / sigma(i)
				diam = (usi .LT. s)

				if (diam) then
					newIndexFound = .true.
					exit
				endif
		
				tied = (usi .EQ. s)
		
				if (tied) then
					mois = (fmois(i) .LT. fm)
					if (mois) then
						newIndexFound = .true.
						exit
					endif
			
					tied = (fmois(i) .EQ. fm)
					if (tied) then
						dens = (dryden (i) .LE. de)
						if (dens) then
							newIndexFound = .true.
							exit
						endif
					endif
				endif
		
				! i is greater than j.
				! Move entry i (the entry we are comparing to) down one:
				sigma(i + 1) = sigma(i)
				fmois(i + 1) = fmois(i)
				dryden(i + 1) = dryden(i)
				key(i + 1) = key(i)
			end do

			if (newIndexFound) then
				! If a new location has been identified record entry j at i + 1 (below).
				newIndexFound = .false. ! Reinitialize for the next comparison.
			else
				i = 0 ! Reseting i to 0 will move this entry (j) to postion 1.
			endif

			! Record the values for fuel index j at the identified index:
			sigma(i + 1) = 1.0 / s
			fmois(i + 1) = fm
			dryden(i + 1) = de
			key(i + 1) = keep
		end do

	end subroutine SORTER

! -- Pagebreak --
! Pg. 92:  This page did not have OCR applied.  I extracted the page and ran OCR on it myself.


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
	subroutine OVLAPS(dryld, sigma, dryden, ak, number, &
						maxno, maxkl, &
						fmois, &
						beta, elam, alone, area)
		implicit none

		! Arguments:
		real*4, intent(in) :: dryld(maxno)			! Ovendry mass per unit area of each element (kg/sq m) (= wdry, ...)
		real*4, intent(in) :: sigma(maxno)			! Surface to volume ratio, 1 / m
		real*4, intent(in) :: dryden(maxno)			! Ovendry mass density, kg / cu m (elsewhere dendry)
		real*4, intent(in) :: ak					! Area influence factor (ak / K_a parameter)
		integer, intent(in) :: number				! The actual number of fuel classes.
		integer, intent(in) :: maxno				! The maximum number of fuel classes allowed.
		integer, intent(in) :: maxkl				! Max triangular matrix size.
		
		real*4, intent(in) :: fmois(maxno)			! Moisture fraction of component
		
		!real*4, intent(out) :: beta(maxno)			! Consolidated interaction matrix (elsewhere = xmat). ! JMR: not quite right!!!!!
		real*4, intent(out) :: beta(maxkl)
		real*4, intent(out) :: elam(maxno, maxno)	! Interaction matrix
		real*4, intent(out) :: alone(maxno)			! Non-interacting fraction for each fuel class.
		real*4, intent(inout) :: area(maxno)		! Fraction of site area expected to be covered at
													! least once by initial planform area of ea size
													! JMR: Change to out!!!!!

		! Locals:
		real :: pi ! Convert to a constant?
		integer :: j, k, l, kj, kl ! Counters
		real :: siga
		real :: a
		real :: bb
		real :: frac ! Intermediate for calculating alone().
		
		!real :: akDyn ! Calculated value of K_a (ak) parameter.
		!real :: akVal ! Value of K_a (ak) parameter, fixed or calculated depending on the mode.
		real :: K_a ! Value of K_a (ak) parameter, fixed or calculated depending on the mode.

		pi = abs(acos(-1.0)) ! Calculate pi.

		! Initialize arrays to 0:
		do j = 1, number
			alone(j) = 0.0
			do k = 1, j
				kj = Loc(j, k)
				beta(kj) = 0.0
			end do
			do k = 1, number
				elam(j, k) = 0.0
			end do
		end do

		do k = 1, number
			!siga = ak * sigma(k) / pi
			do l = 1, k
				
				!kl = Loc(k, l)
				
				! Calculate ak from the fuel moisture of the smaller or similar fuel member (l):
				! Albini & Reinhardt 1997 Equation 4: K_a = K exp(-B * M^2)
				!akDyn = 3.25 * exp(-20 * fmois(l)**2)
				
				! SAV / pi = diameter (units are carried by ak)
				!siga = ak * sigma(k) / pi
				!siga = akDyn * sigma(k) / pi
				
				if (ak .gt. 0) then! or .ge. ?
					! If a valid value has been supplied use a fixed K_a as in the original Burnup:
					K_a = ak
				else
					! A negative value indicates that the K_a should be calculated:
					! Calculate ak from the fuel moisture of the smaller or similar fuel member (l):
					! Albini & Reinhardt 1997 Equation 4: K_a = K exp(-B * M^2)
					K_a = 3.25 * exp(-20 * fmois(l)**2)
				end if
				
				siga = K_a * sigma(k) / pi
				
				kl = Loc(k, l)
				a = siga * dryld(l) / dryden(l) ! siga * ? units in meters
				if (k .EQ. 1) then
					bb = 1.0 - exp(-a)			! JMR: FOFEM suggests this can hit 0?
					area (k) = bb				! JMR: Transcription!!!!
				else
					bb = min(1.0, a)
				end if
				beta(kl) = bb
			end do
		end do

		if (number .EQ. 1) then
			elam (1, 1) = beta (2)				! JMR: Transcription!!!!
			alone(1) = 1.0 - elam(1, 1)
			return
		end if

		do k = 1, number


! -- Pagebreak --
! Pg. 93:


			frac = 0.0
			do l = 1, k
				kl = Loc(k, l)
				frac = frac + beta(kl)
			end do
			if (frac .GT. 1.0) then
				do l = 1, k
					kl = Loc(k, l)
					elam(k, l) = beta(kl) / frac
				end do
				alone(k) = 0.0
			else
				do l = 1, k
					kl = Loc(k, l)
					elam(k, l) = beta(kl)
				end do
				alone(k) = 1.0 - frac
			end if
		end do

		! JMR: Temporary reporting:
		! print *, "OVLAPS()"
! 		print *, "beta:"
! 		print *, beta
! 		print *, "alone:"
! 		print *, alone
! 		print *, "elam:"
! 		print *, elam

	end subroutine OVLAPS


! -- Pagebreak --
! Pg. 94:


	!c This routine initializes variables prior to starting sequence of calls
	!c to subroutine STEP.  On input here, fi is area intensity of spreading
	!c fire, dt is the residence time for the spreading fire.  Only physical
	!c parameters specified are the fuel array descriptors. To call STEP,
	!c one must initialize the following variables.
	! 
	! History: Modernized original Burnup subroutine.
	subroutine START(dt, mxstep, now, maxno, number, wo, alfa, &
						dendry, fmois, cheat, condry, diam, tpig, &
						tchar, xmat, tpamb, fi, flit, fout, &
						tdry, tign, tout, qcum, tcum, acum, qdot, &
						ddot, wodot, work, u, d, r0, dr, &
						ncalls, maxkl)
		implicit none


		! Arguments:
		! JMR_NOTE: The original comments imply that alfa, diam, and wo should all be intent(in).
		! However the code is not consistant with that.
		real*4, intent(in) :: dt				! Spreading fire residence time (s) (= ti, tis, or time elsewhere).
		integer, intent(in) :: mxstep			! Max dimension of historical sequences
		integer, intent(in) :: now 				! Index marks end of time step
		integer, intent(in) :: maxno			! Max number of fuel components
		integer, intent(in) :: number			! Actual number of fuel components
		real*4, intent(inout) :: wo(maxkl)		! Current ovendry loading for the larger of
												! each component pair, kg / sq m.  Updated on return.
		real*4, intent(out) :: alfa(maxno)		! Dry thermal diffusivity of component, sq m / s
		real*4, intent(in) :: dendry(maxno)		! Ovendry density of component, kg /cu m
		real*4, intent(in) :: fmois(maxno)		! Moisture fraction of component
		real*4, intent(in) :: cheat(maxno)		! Specific heat capacity of component, J / kg K
		real*4, intent(in) :: condry(maxno)		! Ovendry thermal conductivity, W / sq m K
		real*4, intent(inout) :: diam(maxkl)	! Current diameter of the larger of each
												! fuel component pair, m.  Updated on return.
		real*4, intent(in) :: tpig(maxno)		! Ignition temperature (K), by component
		real*4, intent(in) :: tchar(maxno)		! tchar = end - pyrolysis temperature (K), by component
		real*4, intent(in) :: xmat(maxkl)		! Table-of-influence fractions between components
		real*4, intent(in) :: tpamb				! Ambient temperature (K)
		real*4, intent(in) :: fi				! Current fire intensity (site avg), kW / sq m

		integer, intent(in) :: maxkl			! Max triangular matrix size.

		! Parameters updated [input and output]:
		integer, intent(inout) :: ncalls		! Counter of calls to this routine
												! = 0 on first call or reset
												! cumulates after first call
												! JMR_NOTE: This is a strange argument as it is only
												! initialized to zero here and is not used.  It is
												! returned and passed on to STEP().  It is probably
												! initialized here to prevent isses related to
												! persistance should more than one simulation be run in
												! on interactive session.
		real*4, intent(out) :: flit(maxno)		! Fraction of each component currently alight
		real*4, intent(out) :: fout(maxno)		! Fraction of each component currently gone out
		real*4, intent(out) :: tdry(maxkl)		! Time of drying start of the larger of each
												! fuel component pair
		real*4, intent(out) :: tign(maxkl)		! Ignition time for the larger of each
												! fuel component pair
		real*4, intent(out) :: tout(maxkl)		! Burnout time of larger component of pairs
		real*4, intent(out) :: qcum(maxkl)		! Cumulative heat input to larger of pair, J / sq m
		real*4, intent(out) :: tcum(maxkl)		! Cumulative temp integral for qcum (drying)
		real*4, intent(out) :: acum(maxkl)		! Heat pulse area for historical rate averaging
		real*4, intent(out) :: qdot(maxkl, mxstep)	! History (post ignite) of heat transfer rate
													! to the larger of each component pair
		real*4, intent(out) :: ddot(maxkl)		! Diameter reduction rate, larger of pair, m / s
		real*4, intent(out) :: wodot(maxkl)		! Dry loading loss rate for larger of pair
		real*4, intent(inout) :: work(maxno)	! Workspace array


! -- Pagebreak --
! Pg. 95:


		! Constant parameters:
		real*4, intent(in) :: u			! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4, intent(in) :: d 		! Fuelbed depth
		real*4, intent(in) :: r0		! Minimum value of mixing parameter
		real*4, intent(in) :: dr		! Max - min value of mixing parameter

		! The constants ch2o and tpdry were included as arguments in the original code.  They have
		! chnaged to globals.
		! The original comments include hvap as a constant, but is not actually used:
		! hvap = heat of vaporization of water J / kg

		! Local constants:
		real, parameter :: rindef = 1.0e+30 ! JMR_NOTE: Make this global?

		! Locals:
		integer :: k, l, kl ! Counters
		real :: delm		! Moisture effect on burning rate (scale factor)
		real :: heatk		! Burn rate factor
		real :: r, tf, ts, thd, tx
		real :: dia			! Diameter for single element.
		real :: cpwet, fac
		real :: dryt		! Drying time for single element.
		real :: tsd
		real :: c			! Thermal conductivity for single element.
		real :: tigk		! Ignition temperature for single element.
		real :: en, e		! Modified Nusselt number obtained from HEATX() (in different places).
		real :: trt			! Minimum ignition time across all fuels (initial estimate?).
		real :: nlit		! Number of elements lit.
		real :: factor
		real :: hb			! "Effective" film heat transfer coefficient returned from HEATX().
		real :: hf			! Film heat transfer coefficient returned from HEATX().
		real :: dtign		! Time to piloted ignition returned from TIGNIT().
		real :: conwet
		real :: aint
		real :: ddt			! Timestep to calculate.  May be less that dt if fuel burns out sooner.
		real :: dnext		! Diameter after combustion in this timestep.
		real :: wnext		! wo after combustion in this timestep.
		real :: df			!

		!c Initialize time varying quantities and, set up work(k)
		!c The diameter reduction rate of fuel component k is given
		!c by the product of the rate of heat tranefer to it, per
		!c unit surface area, and the quantity work(k)

		do k = 1, number
			fout(k) = 0.0
			flit(k) = 0.0
			alfa(k) = condry(k) / (dendry(k) * cheat(k))
			!c effect of moisture content on burning rate (scale factor)
			delm = 1.67 * fmois(k)
			!c effect of component mass density (empirical)
			heatk = dendry(k) / 446.0
			!c empirical burn rate factor, J / cu m - K
			heatk = heatk * 2.01e+06 * (1.0 + delm)
			!c normalize out driving temperature difference (Tfire - Tchar)
			!c to average value of lab experiments used to find above constants
			work(k) = 1.0 / (255.0 * heatk)
			do l = 0, k
				kl = Loc(k, l)
				tout(kl) = rindef			! JMR_NOTE: These are initialized to blank values so they are not inout??????!!!!!
				tign(kl) = rindef
				tdry(kl) = rindef
				tcum(kl) = 0.0
				qcum(kl) = 0.0
			end do
		end do

		! The original code does not initialize acum and qdot.  Failure to do so leads to small
		! variations between runs and changes to results.
		acum = 0.0
		qdot = 0.0
		! There are a large number of other locals that are not explictly initialized.  Testing was
		! done to confrim that explicit initialization was not needed.


! -- Pagebreak --
! Pg. 96: This page did not have OCR applied.  I extracted the page and ran OCR on it myself.


		!c Make first estimate of drying start times for all components
		!c These times are usually brief and make little or no difference

		r = r0 + 0.25 * dr
		tf = tempf(fi, r, tpamb)
		ts = tpamb
		if (tf .LE. (tpdry + 10.)) stop' Igniting fire cannot dry fuel'
		thd = (tpdry - ts) / (tf - ts)
		tx = 0.5 * (ts + tpdry)

		do k = 1, number
			factor = dendry(k) * fmois(k)
			conwet = condry(k) + 4.27e-04 * factor
			do l = 0, k
				kl = Loc(k, l)
				dia = diam (kl)				! JMR: Whitespace!!!!!
				call heatx(u, d, dia, tf, tx, hf, hb, conwet, en)
				call DRYTIM(en, thd, dryt)
				cpwet = cheat(k) + fmois(k) * ch2o
				fac = ((0.5 * dia) ** 2) / conwet
				fac = fac * dendry(k) * cpwet
				dryt = fac * dryt
				tdry(kl) = dryt
				
				! JMR_TEMP_REPORTING:
				if (kl .eq. klVal) then ! Loc(7,1) = 29
					print *, "k, l =", k, l
					print *, "tdry(kl) assignment 1: START()"
					print *, "tdry(kl)", tdry(kl)
					print *, "fac", fac
				end if ! JMR_TEMP_REPORTING end
				
			end do
		end do

		!c Next, determine which components are alight in spreading fire

		tsd = tpdry

		do k = 1, number
			c = condry(k)
			tigk = tpig(k)
			do l = 0, k
				kl = Loc(k, l)
				dryt = tdry(kl)
				if (dryt .LT. dt) then
					dia = diam (kl)									! JMR: Whitespace!!!!!
					ts = 0.5 * (tsd + tigk)
					call heatx(u, d, dia, tf, ts, hf, hb, c, e)
					tcum(kl) = max((tf - ts) * (dt - dryt), 0.)		! JMR: Trailing 0s!!!!!
					qcum(kl) = hb * tcum(kl)
					if (tf .GT. (tigk + 10.)) then
						call TIGNIT(tpamb, tpdry, tpig(k), tf, condry(k), &
									 cheat(k), fmois(k), dendry(k), hb, dtign)
						trt = dryt + dtign
						tign(kl) = 0.5 * trt
						
						! JMR_TEMP_REPORTING:
						!if (k .eq. 7) then ! Fuel type 7 is the first fuel showing issues.
						!if (kl .eq. 29) then! Loc(7,1) = 29
						if (kl .eq. klVal) then
							print *, "tign(kl) assignment 1: START()"
							! Time should be the first timestep.
							print *, "tign(kl)", tign(kl)
							!print *, "tpamb", tpamb
							!print *, "tpdry", tpdry
							!print *, "tpig(k)", tpig(k)
							print *, "tf", tf
							!print *, "condry(k)", condry(k)
							!print *, "cheat(k)", cheat(k)
							!print *, "fmois(k)", fmois(k)
							!print *, "dendry(k)", dendry(k)
							print *, "hb", hb
							print *, "dtign", dtign
							print *, "dryt", dryt
						end if ! JMR_TEMP_REPORTING end
						
						if (dt .GT. trt) then
							flit(k) = flit(k) + xmat(kl)
						end if
					end if
				end if
			end do
		end do

		nlit = 0


! -- Pagebreak --
! Pg. 97:


		trt = rindef

		!c Determine minimum ignition time and verify ignition exists

		do k = 1, number
			if (flit(k) .GT. 0.) nlit = nlit + 1			! JMR: Oneliner!!!!!
			do l = 0, k
				kl = Loc(k, l)
				trt = min(trt, tign(kl))
			end do
		end do

			if (nlit .EQ. 0) stop' START ignites no fuel'		! JMR: Oneliner!!!!!, indenting!

		!c Deduct trt from all time estimates, resetting time origin

		do k = 1, number
			do l = 0, k
				kl = Loc(k, l)
				if (tdry(kl) .LT. rindef) then
					tdry(kl) = tdry (kl) - trt				! JMR: Whitespace!!!!!
					
					! JMR_TEMP_REPORTING:
					!if (kl .eq. 29) then ! Loc(7,1) = 29
					if (kl .eq. klVal) then
						print *, "tdry(kl) assignment 2: START()"
						print *, "tdry(kl)", tdry(kl)
						print *, "trt", trt
					end if ! JMR_TEMP_REPORTING end
					
				end if
				if (tign(kl) .LT. rindef) then
					tign(kl) = tign(kl) - trt
				end if
			end do
		end do

		!c Now go through all component pairs and establish burning rates
		!c for all the components that are ignited; extrapolate to end time dt

		do k = 1, number
			if (flit(k) .EQ. 0.) then
				do l = 0, k
					kl = Loc(k, l)
					ddot(kl) = 0.0
					tout(kl) = rindef
					wodot(kl) = 0.0
				end do
			else
				ts = tchar(k)
				c = condry(k)
				do l = 0, k
					kl = Loc(k, l)
					dia = diam(kl)
					call heatx(u, d, dia, tf, ts, hf, hb, c, e)
					qdot(kl, now) = hb * max((tf - ts), 0.)			! JMR: Trailing 0s!!!!!
					aint = (c / hb) ** 2
					ddt = dt - tign(kl)
					acum(kl) = aint * ddt
					ddot(kl) = qdot(kl, now) * work(k)
					tout(kl) = dia / ddot(kl)
					dnext = max(0., (dia - ddt * ddot(kl)))			! JMR: Trailing 0s!!!!!
					wnext = wo(kl) * ((dnext / dia) ** 2)
					wodot(kl) = (wo(kl) - wnext) / ddt


! -- Pagebreak --
! Pg. 98:


					diam(kl) = dnext
					wo(kl) = wnext
					df = 0.0
					if (dnext .LE. 0.) then			! JMR: Trailing 0s!!!!!
						df = xmat(kl)
						wodot(kl) = 0.0
						ddot(kl) = 0.0
					end if
					flit(k) = flit(k) - df
					fout(k) = fout(k) + df
				end do
			end if
		end do

		ncalls = 0

	end subroutine START


! -- Pagebreak --
! Pg. 99:


	!c Computes fi = site avg fire intensity given the burning rates of all
	!c interacting pairs of fuel components [ wodot ], the mineral ash content
	!c of each component [ ash ], the heat of combustion value [ htval ] for
	!c each, and the number of fuel components [ number ], where max - maxno.
	!c fi is in kW / sq m, while htval is in J / kg.
	!
	!c fint(k) is the correction to fi to adjust
	!c the intensity level to be the local value where size k is burning.
	!
	! History: Modernized original Burnup subroutine.
	subroutine FIRINT(wodot, ash, htval, maxno, number, maxkl, &
						area, fint, fi)
		implicit none

		! Arguments:
		real*4, intent(in) :: wodot(maxkl)	! Burning rates of interacting pairs of fuel components
		real*4, intent(in) :: ash(maxno)		! Mineral content, fraction dry mass
		real*4, intent(in) :: htval(maxno)	! Low heat of combustion, J / kg
		integer, intent(in) :: maxno			! The maximum number of fuel classes allowed.
		integer, intent(in) :: number			! The actual number of fuel classes.
		integer, intent(in) :: maxkl			! Max triangular matrix size.
		real*4, intent(in) :: area(maxno)		! Fraction of site area expected to be covered at
												! least once by initial planform area of ea size
		real*4, intent(out) :: fint(maxno)	! Corrected local fire intensity for each fuel type.
		real*4, intent(out) :: fi				! Site avg fire intensity (kW / sq m)

		! Locals:
		real :: sum ! Running total for fi.
		integer :: k, l, kl ! Counters
		real :: wdotk
		real :: term
		real :: ark ! Area for element k.

		! Constants:
		real, parameter :: small = 1.0e-06

		sum = 0.0
		do k = 1, number
			wdotk = 0.0
			do l = 0, k
				kl = Loc(k, l)
				wdotk = wdotk + wodot(kl)
			end do
			term = (1.0 - ash(k)) * htval(k) * wdotk * 1.0e-03
			ark = area(k)
			if (ark .GT. small) then
				fint(k) = term / ark - term
			else
				fint(k) = 0.0
			end if
			sum = sum + term
		end do

		fi = sum

	end subroutine FIRINT

! -- Pagebreak --
! Pg. 100:



	!c This routine stashes output from the BURNUP model package on a snapshot
	!c basis.  Every time it is called, it "dumps" a picture of the status of
	!c each component of the fuel complex, as a table of interacting pairs.
	! This routine leaves the history file open when it returns.  It does not return it's unit
	! identifier (mum) so the main program must know that number, which is not very robust.
	!
	! History: Modernized original Burnup subroutine.
	subroutine STASH(time, now, maxno, number, outfil, fi, &
						flit, fout, wo, wodot, diam, ddot, &
						tdry, tign, tout, fmois, maxkl, nun)
		implicit none

		! Arguments:
		real, intent(in) :: time 				! Igniting surface fire residence time (s)
		integer, intent(in) :: now				! Time step?
		integer, intent(in) :: maxno			! The maximum number of fuel classes allowed.
		integer, intent(in) :: number			! The actual number of fuel classes.
		character*12, intent(out) :: outfil		! The name of the output file.
		! NOTE: There seems to be no reason why outfil is returned.  It appears to be overwritten
		! wherever it is used subsequent to this call and could be safely changed to a local
		! variable.  It is maintained only to preserve the original interface.
		real*4, intent(in) :: fi				! Site avg fire intensity (kW / sq m)
		real*4, intent(in) :: flit(maxno)		! Fraction of each component currently alight
		real*4, intent(in) :: fout(maxno)		! Fraction of each component currently gone out
		real*4, intent(in) :: wo(maxkl)			! Current ovendry loading for the larger of
												! each c-omponent pair, kg / sq m
		real*4, intent(in) :: wodot(maxkl)		! Burning rates of interacting pairs of fuel components
		real*4, intent(in) :: diam(maxkl)		! Current diameter of the larger of each fuel component pair (m)
		real*4, intent(in) :: ddot(maxkl)		! Diameter reduction rate, larger of pair, m / s
		real*4, intent(in) :: tdry(maxkl)		! Time of drying start of the larger of each fuel component pair
		real*4, intent(in) :: tign(maxkl)		! Ignition time for the larger of each fuel component pair
		real*4, intent(in) :: tout(maxkl)		! Burnout time of larger component of pairs
		real*4, intent(in) :: fmois(maxno)		! Moisture fraction of component
		integer, intent(in) :: maxkl			! Max triangular matrix size.
		integer, intent(inout) :: nun			! Stash file unit identifier, returned on first call,
												! passed in for subsequent.

		! Locals:
		!character*12 :: outfil
		character*12 :: histry
		logical snaps

		integer :: m, n, mn ! Counter.
		real :: wd, wg, wdm, wgm ! Accumulators
		real :: wd0, wg0													! JMR: Are these are persistant?????
		real :: wgf, wdf ! Ratios?
		real :: fmm ! Individual fuel moisture
		integer :: istash ! User return value: Should stats be output (stashed) at each timepoint.
		integer :: mum ! History file unit identifier.
		integer :: ido ! User menu return value.
		integer :: openStat, writeStat ! IO status.

		! Constants:
		character (len = *), parameter :: outFormat = "(8e10.3)"

		! If this is the first time step:
		if (now .EQ. 1) then
			! Calculate wd and wg at time 0:
			! JMR_NOTE: These will need to be preserved in C++ if persistant?????
			wd = 0.0
			wg = 0.0
			do m = 1, number
				fmm = fmois(m)
				wdm = 0.0
				do n= 0, m					! JMR: Whitespace
					mn = Loc(m, n)
					wdm = wdm + wo(mn)
				end do
				wgm = wdm * (1.0 + fmm)
				wd = wd + wdm
				wg = wg + wgm
			end do
			wd0 = wd
			wg0 = wg
	
			nun = 77
	
			! Open the history file:
			mum = 66
			histry = 'HISTORY.DAT'
			open(mum, file=histry, status='UNKNOWN', form='FORMATTED')
	
			! Determine if the user if they want to record snapshots:
			snaps = .FALSE.
			do ! Loop until valid input is received. Alt: Set istash = -1, do while (istash .ne. 0 .and. istash .ne. 1)
				write(*, "(' Enter 1 to stash fuelbed snapshots, O to skip ', $)")
				read(*, *) istash
				if (istash .EQ. 0) then
					exit
				else if (istash .EQ. 1) then
					snaps = .TRUE.
					exit ! Continue to the stash file setup.
				end if
				! Otherwise start again.
			end do
	
			! This could be moved into the above loop but I'm leaving it here to keep as much of the original structure as possible:
			if (istash .EQ. 1) then
				! Prompt the user for a output file name and open it:
				do ! Name loop: Loop until a valid file is specified.
					write(*,"(' File name [ max 12 char ] for STASH output    ', $)")
					read(*, fmt = '(a12)') outfil

					open(nun, file = outfil, status = 'NEW', form = 'FORMATTED', iostat = openStat)
			
					! If there is an error opening the file:
					do while (openStat .ne. 0)
						write(*, "(' Error opening file    'a12)") outfil
						write(*, *) 'Enter 0 to try a new name'
						write(*, *) 'Enter 1 to overwrite existing file'

! -- Pagebreak --
! Pg. 101:

						write(*,  "(' Enter 2 to terminate program now    ', $)")

						read(*,*) ido
						if (ido .EQ. 0) then
							exit

						else if (ido .EQ. 1) then
							open(nun, file = outfil, status = 'NEW', form = 'FORMATTED', iostat = openStat)

						else if (ido .EQ. 2) then
							stop' Program ended'
						!else ! Implied.
						!	cycle
						end if
					end do ! File error loop.
			
					if (openStat .eq. 0) then
						exit
					end if
					! If there is still an error at this point a new name is needed:
				end do ! File name loop.
			end if ! (istash .EQ. 1)
		end if ! (now .EQ. 1)

		if (snaps) then
			write(nun, outFormat) time, fi
		end if

		wd = 0.0
		wg = 0.0
		do m = 1, number
			fmm = fmois(m)
			wdm = 0.0
			if (snaps) then
				write(nun, outFormat) time, flit(m), fout(m)
			end if
			do n = 0, m
				mn = Loc(m, n)
				wdm = wdm + wo(mn)
				if (snaps) then
					write(nun, outFormat) time, wo(mn), &
							wodot(mn), diam(mn), ddot(mn), &
							tdry(mn), tign(mn), tout(mn)
				end if

			end do
			wgm = wdm * (1.0 + fmm)
			wd = wd + wdm
			wg = wg + wgm
		end do
		wgf = wg / wg0
		wdf = wd / wd0

		write(mum, outFormat) time, wg, wd, wgf, wdf, fi

	end subroutine STASH


! -- Pagebreak --
! Pg. 102: This page did not have OCR applied.  I extracted the page and ran OCR on it myself.


	! This routine computes the halfspace surface ignition time under steady radiant heating with
	! surface film cooling.
	! History: Modernized original Burnup subroutine.
	!
	! JMR_NOTE: If this only has one return value it could be turned into a function.
	subroutine TIGNIT(tpam, tpdr, tpig, tpfi, cond, &
						chtd, fmof, dend, hbar, tmig)
		implicit none

		! Arguments:
		real*4, intent(in) :: tpam	! ambient temperature, K
		real*4, intent(in) :: tpdr	! fuel temperature at start of dryirig, K			! JMR: Spelling!!!!!
									! Currently this is always tpdry, so this argument could be cut.
		real*4, intent(in) :: tpig	! fuel surface temperature at ignition, K
		real*4, intent(in) :: tpfi	! fire enviroriment temperature, K
		real*4, intent(in) :: cond	! fuel ovendry thermal conductivity, W / m K
		real*4, intent(in) :: chtd	! fuel ovendry specific heat capacity, J / kg K
		real*4, intent(in) :: fmof	! fuel moisture content, fraction dry weight
		real*4, intent(in) :: dend	! fuel ovendry density, kg / cu m
		real*4, intent(in) :: hbar	! effective film heat transfer coefficient [< HEATX] W / sq m K
		real*4, intent(out) :: tmig	! predicted time to piloted ignition, s

		! Locals:
		real :: b03
		real :: xlo, xhi, xav	! Binary search bounds and middle (average).
		real :: fav				! Value of ff() for current search value.
		real :: beta, conw
		real :: dtb				! The temperature increase required to reach the drying temperature.
		real :: dti				! The temperature increase required to reach the ignition temperature.
		real :: ratio, rhoc

		! Constants:
		real, parameter :: a03 = -1.3371565	! ff() parameter 1
		real, parameter :: a13 = 0.4653628	! ff() parameter 2
		real, parameter :: a23 = -0.1282064	! ff() parameter 3
		real, parameter :: pinv = 2.125534
		real, parameter :: small = 1.e-06	! Terminate the search when we are this close.
		real, parameter :: hvap = 2.177e+06 ! Heat of vaporization of water J/kg.
		real, parameter :: cpm = 4186.0
		real, parameter :: conc = 4.27e-04

		!c radiant heating equivalent form gives end condition fixes beta

		b03 = a03 * (tpfi - tpig) / (tpfi - tpam)

		!c find x that solves ff(x) = 0 ; method is binary search

		xlo = 0.0
		xhi = 1.0
		! Testing was done to confrim that explicit initialization of other locals was not needed.

		do
			xav = 0.5 * (xlo + xhi)
	
			! The original code implements this as a statement function:
			!fav = ff(xav)
			! Original notes:
			!c approximate function of beta to be solved is ff(x) where
			!c x  =  1 / (1 + p * beta)  { Hastings, Approximations for
			!c digital computers } and we employ      pinv  =  1 / p
			fav = b03 + xav * (a13 + xav * (a23 + xav))
			! Or:
			!fav = b03 + xav * (a13 + (a23 * xav + xav * xav))
			!fav = b03 + (a13 * xav) + (a23 * xav ** 2) + (xav ** 3)
	
			if (abs(fav) .LE. small) then
				exit
			else if (fav .LT. 0.) then		! JMR: Trailing 0s!!!!!
				xlo = xav
			else if (fav .GT. 0.) then
				xhi = xav
			end if
		end do

		beta = pinv * (1.0 - xav) / xav
		conw = cond + conc * dend * fmof
		dtb = tpdr - tpam
		dti = tpig - tpam
		ratio = (hvap + cpm * dtb) / (chtd * dti)
		rhoc = dend * chtd * (1.0 + fmof * ratio)
		tmig = ((beta / hbar) ** 2) * conw * rhoc

	end subroutine TIGNIT


! -- Pagebreak --
! Pg. 103:


	!c Given a Nusselt number (enu, actually Biot number = h D / k)
	!c and the dimensionless temperature rise required for the start
	!c of surface drying (theta), returns the dimerisionless time (tau)
	!c needed to achieve it. The time is given multiplied by thermal
	!c diffusivity and divided by radius squared. Solution by binary search.
	!
	! History: Modernized original Burnup subroutine.
	subroutine DRYTIM(enu, theta, tau)
		implicit none

		! Arguments:
		real, intent(in) :: enu		! Biot number.
		real, intent(in) :: theta	! Temperature rise required for the start of moisture loss.
		real, intent(out) :: tau	! Time required for the start of moisture loss.

		! Locals:
		real :: xl, xh, xm ! The binary search low and high bounds, and the search center.
		real :: x
		real*4 :: approx
		integer n ! Counter

		! Constants:
		real, parameter :: p = 0.47047

		xl = 0.0
		xh = 1.0

		! No rationale is given for the use of 15 cycles in the binary search.  I assume that is was
		! determined to be sufficient empirically.  Using a limit for convergence might be better.
		! Additionally, since the search middle is calculated at the start ot the loop only 14
		! cycles actually inform the result.
		do n = 1, 15
			xm = 0.5 * (xl + xh)
			approx = ErrorApprox(xm, theta)
			if (approx .LT. 0.) then ! or if (ErrorAppox(xm) .LT. 0.) then
				xl = xm
			else
				xh = xm
			end if
		end do

		x = (1.0 / xm - 1.0) / p
		tau = (0.5 * x / enu) **2

	end subroutine DRYTIM


! -- Pagebreak --
! Pg. 104:


	! This routine calculates how heat is transfered from the fire environment to a given fuel type.
	!c Given horizontal windspeed u at height d [top of fuelbed], cylindrical
	!c fuel particle diameter dia, fire environment temperature tf, and mean
	!c surface temperature, ts, subroutine returns film heat transfer coefficient
	!c hfm and an "effective" film heat transfer coefficient including radiation
	!c heat transfer, hbar.  Using the wood's thermal conductivity, cond, the
	!c modified Nusselt number [ en ] used to estimate onset of surface drying
	!c is returned as well.
	!
	! History: Modernized original Burnup subroutine.
	subroutine HEATX(u, d, dia, tf, ts, hfm, hbar, cond, en)
		implicit none

		! Arguments:
		real*4, intent(in) :: u		! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4, intent(in) :: d		! Fuelbed depth (m)
		real*4, intent(in) :: dia	! Fuel diameter
		real*4, intent(in) :: tf	! Fire environment temperature
		real*4, intent(in) :: ts	! Mean surface temperature
		real*4, intent(out) :: hfm	! Film heat transfer coefficient
		real*4, intent(out) :: hbar	! "Effective" film heat transfer coefficient
		real*4, intent(in) :: cond	! Wood thermal conductivity
		real*4, intent(out) :: en	! Modified Nusselt number

		! Locals:
		real :: v		! Estimate of relative vertical air velocity over fuel element 
		real :: re		! Reynolds number (air)
		real :: enuair	! Nusselt number
		real :: conair	! (Forced) convection of air?
		real :: fac
		real :: hfmin	! Film heat transfer coefficient for natural convection (used as minimum value)
		real :: hrad	! Radiation contribution

		! Constants:
		real, parameter :: g = 9.8
		real, parameter :: vis = 7.5e-05	! Kinematic viscosity of hot air
		real, parameter :: a = 8.75e-03
		real, parameter :: b = 5.75e-05
		real, parameter :: rad = 5.67e-08	! Stefan-Boltzmann radiation constant(W/m^2-K^4)
		real, parameter :: fmfac = 0.382
		real, parameter :: hradf = 0.5		! View factor emissivity

		hfm = 0.0
		! Testing was done to confrim that explicit initialization of other locals was not needed.

		if (dia .gt. b) then
			v = sqrt(u * u + 0.53 * g * d)
			re = v * dia / vis
			enuair = 0.344 * (re**0.56)
			conair = a + b * tf
			fac = sqrt(abs(tf - ts) / dia)
			hfmin = fmfac * sqrt(fac)
			hfm = max((enuair * conair / dia), hfmin)
		endif

		hrad = hradf * rad * (tf + ts) * (tf * tf + ts * ts)
		hbar = hfm + hrad
		en = hbar * dia / cond

	end subroutine HEATX


! -- Pagebreak --
! Pg. 105:


	!c Returns a fire environment temperature, TEMPF, given the fire intensity
	!c q in kW / square meter, the ambient temperature tamb in Kelvins, and the
	!c dimensionless mixing parameter r.
	!
	! History: Modernized original Burnup subroutine.
	function TEMPF(q, r, tamb)
		implicit none

		! Arguments:
		real, intent(in) :: q		! Fire intensity
		real, intent(in) :: r		! Dimensionless mixing parameter
		real, intent(in) :: tamb	! Ambient temperature (K)

		! Locals:
		real :: term, rlast, den, rnext
		real*4 :: tempf !TEMPF ! Return value.

		! Constants:
		real, parameter :: err = 1.0e-04
		real, parameter :: aa = 20.0

		! Testing was done to confrim that explicit initialization of locals was not needed here.

		term = r / (aa * q)
		rlast = r

		do
			den = 1.0 + term * (rlast + 1.0) * (rlast * rlast + 1.0)
			rnext = 0.5 * (rlast + 1.0 + r / den)
			if (abs(rnext - rlast) .LT. err) then
				tempf = rnext * tamb
				return
			end if
			rlast = rnext
		end do

	end function TEMPF


! -- Pagebreak --
! Pg. 106:


	! This routine calculates one timestep of the fuel consumption process.
	!
	! History: Modernized original Burnup subroutine.
	! JMR_NOTE: This routine takes a large number of arguments and the order is a bit confusing
	! with input and output parameters mixed in the order.
	subroutine STEP(dt, mxstep, now, maxno, number, wo, alfa, &
					dendry, fmois, cheat, condry, diam, tpig, &
					tchar, xmat, tpamb, fi, flit, fout, &
					tdry, tign, tout, qcum, tcum, acum, qdot, &
					ddot, wodot, work, u, d, r0, dr, &
					ncalls, maxkl, tin, fint, fid)
		implicit none

		!c Updates status of all fuel component pairs and returns a snapshot

		!c Input parameters:

! JMR_NOTE: All explicitly declared reals were real*4 but can probably be changed....
		! JMR_NOTE: These are in original order, not argument order.
		real*4, intent(in) :: dt				! time step, sec
		integer, intent(in) :: mxstep			! max dimension of historical seouences
		integer, intent(in) :: now				! index marks end of time step
		integer, intent(in) :: maxno			! max number of fuel components
		integer, intent(in) :: number			! actual number of fuel components
		real*4, intent(inout) :: wo(maxkl)		! current ovendry loading for the larger of	
												! each component pair, kg / sq m
		real*4, intent(in) :: alfa(maxno)		! dry thermal diffusivity of component, sq m / s
		real*4, intent(in) :: dendry(maxno)		! ovendry density of component, kg / cu m
		real*4, intent(in) :: fmois(maxno)		! moisture fraction of component
		real*4, intent(in) :: cheat(maxno)		! specific heat capacity of component, J / kg K
		real*4, intent(in) :: condry(maxno)		! ovendry thermal conductivity, w / sq m K
		real*4, intent(inout) :: diam(maxkl)	! current diameter of the larger of each			! JMR: Intent?????
												! fuel component pair, m
		real*4, intent(in) :: tpig(maxno)		! ignition temperature (K), by component
		real*4, intent(in) :: tchar(maxno)		! end - pyrolysis temperature (K), by component
		real*4, intent(in) :: xmat(maxkl)		! table of influence fractions between components
		real*4, intent(in) :: tpamb				! ambient temperature (K)
		real*4, intent(in) :: fi				! current fire intensity (site avg), kW / sq m
		real*4, intent(in) :: work(maxno)		! factor of heat transfer rate hbar * (Tfire - Tebar)
												! that yields ddot (k)
												! plus the following constants and bookk keeping parameters
												! u, d, r0, dr, ch20, ncalls, maxkl

		! Not in argument order...
		!c Constant parameters
		real*4, intent(in) :: u 				! Mean horizontal windspeed at top of fuelbed (m/s).
		real*4, intent(in) :: d					! Fuelbed depth (m)
		real*4, intent(in) :: r0				! minimum value of mixing parameter
		real*4, intent(in) :: dr				! max - min value of mixing parameter

		integer, intent(in) :: maxkl			! Max triangular matrix size.
		real*4, intent(in) :: tin				! start of current time step
		real*4, intent(in) :: fint(maxno)		! correction to fi to compute local intensity
												! that may be different due to k burning
		real*4, intent(in) :: fid				! fire intensity due to duff burning ... this is
												! used to up the fire intensity for fuel pieces
												! that are burning without interacting with others


		!c Parameters updated [input and output)

		integer, intent(inout) :: ncalls		! counter of calls to this routine ... 
												! = 0 on first call or reset
												! cumulates after first call
		real*4, intent(inout) :: flit(maxno)	! fraction of each component currently alight
		real*4, intent(inout) :: fout(maxno)	! fraction of each component currently gone out
		real*4, intent(inout) :: tdry(maxkl)	! time of drying start of the larger of each
												! fuel component pair
		real*4, intent(inout) :: tign(maxkl)	! ignition time for the larger of each
												! fuel component pair
		real*4, intent(inout) :: tout(maxkl)	! burnout time of larger component of pairs
		real*4, intent(inout) :: qcum(maxkl)	! cumulative heat input to larger of pair, J / sq m

		! The constants ch2o and tpdry were included as arguments in the original code.  They have
		! chnaged to globals.
		! Note: The original code documents the following variable, but it is mot actually used.
		!real*4, intent(in) :: hvap			! heat of vaporization of water, J / kg


! -- Pagebreak --
! Pg. 107:


		real*4, intent(inout) :: tcum(maxkl)	! cumulative temp integral for qcum (drying)
		real*4, intent(inout) :: acum(maxkl)	! heat pulse area for historical rate averaging
		real*4, intent(inout) :: qdot(maxkl, mxstep)	! history (post ignite) of heat transfer rate
													! to the larger of component pair, W / sq m
		real*4, intent(inout) :: ddot(maxkl)	! diameter reduction rate, larger of pair, m / s
		real*4, intent(inout) :: wodot(maxkl)	! dry loading loss rate for larger of pair

		! Locals: (not in a consistant order.)
		logical flag
		real :: tnow, tnext ! The time of this and the next timestep.
		real :: tdun	! The burnout time for a single pair.
		real :: tgo		! Time left to burnout.
		real :: tifi	! Time when fire ignition phase ended
		real :: next
		real :: gi
		integer :: nspan
		real :: tst, aint, qqq
		real :: tav1, tav2, tav3, tavg ! Time over which to perform averaging
		real :: tbar
		integer :: index
		real :: qdsum	! Sum of heat transfer (W/m^2 * s = J/m^2 ?).
		real :: qdavg	! Average heat transfer...
		real :: deltim, rate, dryt, dqdt
		real :: qd
		real :: dteff, heff, delt, factor, dtef
		real :: he		! qcum / tcum
		real*4 :: tf, ts
		real :: biot
		real :: cpwet
		real :: c		! Thermal conductivity for a single fuel component.
		real :: conwet
		real :: ddt		! Timestep to calculate.  May be less that dt if fuel burns out sooner.
		real :: dia		! Diameter for a single fuel component (kl).
		real :: dnext	! Diameter after combustion in this timestep.
		real :: wnext	! wo after combustion in this timestep.
		real :: dtcum
		real :: dtlite	! Time to ignition returned by TIGNIT().
		real :: e
		real :: fac
		real :: hb, hf
		real :: r
		real :: tfe
		real :: thd
		real :: tlit	! Ignition time for a single fuel component.
		real :: tspan
		real :: dtemp

		integer :: k, l, mu, kl ! Counters (kl is a bit different)

		! Constants:
		real, parameter :: rindef = 1.0e+30

		! There are a large number of locals in this routine that are not explictly initialized.
		! Testing was done to confrim that explicit initialization was not needed.

		ncalls = ncalls + 1
		tnow = tin
		tnext = tnow + dt
		!c tifi = time when fire ignition phase ended (at now = 1)
		tifi = tnow - float(now - 1) * dt
		next = now + 1

		kLoop : do k = 1, number
			c = condry(k)
			lLoop : do l = 0, k
				kl = Loc(k, l)
				tdun = tout(kl)

				!c See if k of (k, l) pair burned out

				if (tnow .GE. tdun) then
					ddot(kl) = 0.0
					wodot(kl) = 0.0
					!goto 10 ! JMR: Remove!!!!!
					cycle lLoop
				end if
				if (tnext .GE. tdun) then
					tgo = tdun - tnow
					ddot(kl) = diam(kl) / tgo


! -- Pagebreak --
! Pg. 108:


					wodot(kl) = wo(kl) / tgo
					wo(kl) = 0.0
					diam (kl) = 0.0 ! JMR: Whitespace!!!!!
					cycle lLoop
				end if

				!c k has not yet burned out ... see if k of (k, 1) pair is ignited ! JMR: Transcription!!!!!

				tlit = tign(kl)
				if (tnow .GE. tlit) then
					ts = tchar(k)
					if (l .EQ. 0) then
						r = r0 + 0.5 * dr
						gi = fi + fid
					end if
					if ((l .NE. 0) .AND. (l .NE. k)) then
						r = r0 + 0.5 * (1.0 + flit(l)) * dr
						gi = fi + fint(k) + flit(l) * fint(l)
					end if								! JMR: Link with following if!!!!!!
					if (l .EQ. k) then
						r = r0 + 0.5 * (1.0 + flit(k)) * dr
						gi = fi + flit (k) * fint(k)
					end if
					tf = tempf(gi, r, tpamb)
					dia = diam(kl)
					call heatx(u, d, dia, tf, ts, hf, hb, c, e)
					qqq = hb * max((tf - ts),  0.)
					tst = max(tlit, tifi)
					nspan = max(l, nint((tnext - tst) / dt))
					if (nspan .LE. mxstep) qdot(kl, nspan) = qqq ! JMR: Oneliner!!!!! Connect with next if...
					if (nspan .GT. mxstep) then
						do mu = 2, mxstep
							!qdot(kl, mu - l) = qdot (kl, mu) ! JMR: Transcription: l for 1 substitution!!!!!
							qdot(kl, mu - 1) = qdot(kl, mu)
						end do
						qdot(kl, mxstep) = qqq
					end if
					aint = (c / hb) ** 2
					acum(kl) = acum(kl) + aint * dt
					
					! Time over which to perform averaging:
					tav1 = tnext - tlit ! Time since ignition.
					tav2 = acum(kl) / alfa(k) ! Measure of square of distance heat has penetrated fuel.
					tav3 = ((dia / 4.) ** 2) / alfa (k) ! Measure of time heat takes to reach center of fuel.
					tavg = min(tav1, tav2, tav3)
					
					!index = l + min(nspan, mxstep) ! JMR: Transcription: l for 1 substitution!!!!!
					index = 1 + min(nspan, mxstep)
					qdsum = 0.0
					tspan = 0.0
					deltim = dt
					
					! JMR: Mod -> parallel!!!!!
					! Calculate qdsum (sum of heat transfer (W/m^2 * s = J/m^2)):
					do
						index = index - 1
						if (index .EQ. 1) then
							deltim = tnext - tspan - tlit
						endif
						
						if ((tspan + deltim) .GE. tavg) then
							deltim = tavg - tspan
						end if
						
						qdsum = qdsum + qdot (kl, index) * deltim ! JMR: Whitespace!!!!!
						tspan = tspan + deltim
						
						if ((tspan .LT. tavg) .AND. (index .GT. 1)) then ! Will this pass the first time?  Yes!
							cycle
						else
							exit
						endif
					end do

					qdavg = max (qdsum / tspan, 0.) ! JMR: trailing zeros!!!!!
					ddot(kl) = qdavg * work(k)
					dnext = max(0., dia - dt * ddot(kl)) ! JMR: trailing zeros!!!!!


! -- Pagebreak --
! Pg. 109:


					wnext = wo(kl) * ((dnext / dia) ** 2)
					if ((dnext .EQ. 0.) .AND. (ddot(kl) .GT. 0.)) then
						tout(kl) = tnow + dia / ddot(kl)
					end if
					if ((dnext .GT. 0.) .AND. (dnext .LT. dia)) then
						rate = dia / (dia - dnext)
						tout(kl) = tnow + rate * dt
					end if
					if (qdavg .LE. 20.) tout(kl) = 0.5*(tnow + tnext) ! JMR: Whitespace...
					ddt = min(dt, (tout(kl) - tnow))
					wodot(kl) = (wo(kl) - wnext) / ddt
					diam(kl) = dnext
					wo(kl) = wnext
					cycle lLoop
				end if ! (tnow .GE. tlit) ! JMR: Mod -> parallel!!!!!

				!c See if k of (k, l) has reached outer surface drying stage yet

				dryt = tdry(kl)
				if ((tnow .GE. dryt) .AND. (tnow .LT. tlit)) then
					if (l .EQ. 0) then
						r = r0
						gi = fi + fid
					end if
					if (l .EQ. k) then
						r = r0
						gi = fi
					end if
					if ((l .NE. 0) .AND. (l .NE. k)) then
						r = r0 + 0.5 * flit(l) * dr
						gi = fi + flit (l) * fint(l)
					end if
					tf = tempf(gi, r, tpamb)
					ts = tpamb
					dia = diam (kl) ! JMR: Transcription!!!!!
					call heatx(u, d, dia, tf, ts, hf, hb, c, e)
					dtemp = max(0., (tf - ts))
					dqdt = hb * dtemp
					qcum(kl) = qcum(kl) + dqdt * dt
					tcum(kl) = tcum(kl) + dtemp * dt
					dteff = tcum(kl) / (tnext - dryt)
					heff = qcum(kl) / tcum(kl)
					tfe = ts + dteff
					dtlite = rindef
			
					if (.not. (tfe .LE. (tpig (k) + 10.))) then
						call TIGNIT(tpamb, tpdry, tpig(k), tfe, &
							condry(k), cheat(k), fmois(k), dendry(k), &
							heff, dtlite)
					endif
					tign(kl) = 0.5 * (dryt + dtlite)
					
					! JMR_TEMP_REPORTING:
					!if (k .eq. 7) then ! Fuel type 7 is the first fuel showing issues.
					!if (kl .eq. 29) then ! Loc(7,1) = 29
					!if (kl .eq. 7) then ! Loc(3,1)
					if (kl .eq. klVal) then
						print *, "tign(kl) assignment 2: STEP()"
						print *, "tnow", tnow
						print *, "tign(kl)", tign(kl)
						!print *, "tpamb", tpamb
						!print *, "tpdry", tpdry
						!print *, "tpig(k)", tpig(k)
						print *, "tfe", tfe
						! Inputs to tfe:
						!print *, "ts", ts
						print *, "dteff", dteff
						print *, "tcum(kl)", tcum(kl)
						!print *, "tnext", tnext
						
! 						!print *, "condry(k)", condry(k)
! 						!print *, "cheat(k)", cheat(k)
! 						!print *, "fmois(k)", fmois(k)
! 						!print *, "dendry(k)", dendry(k)
						print *, "heff", heff
						! Inputs to heff:
						! tcum(kl) is above.
						print *, "qcum(kl)", qcum(kl)
						
						print *, "dtlite", dtlite
						print *, "dryt", dryt
					end if ! JMR_TEMP_REPORTING end

					!c If k will ignite before time step over, must interpolate

					if (tnext .GT. tign(kl)) then
						ts = tchar(k)
						call heatx(u, d, dia, tf, ts, hf, hb, c, e)


! -- Pagebreak --
! Pg. 110:


						qdot(kl, 1) = hb * max((tf - ts), 0.)
						qd = qdot(kl, 1)
						ddot(kl) = qd * work(k)
						delt = tnext - tign(kl)
						dnext = max(0., dia - delt * ddot(kl))
						wnext = wo(kl) * ((dnext / dia) ** 2)
						if (dnext .EQ. 0.) tout(kl) = tnow + dia / ddot(kl)
						if ((dnext .GT. 0.) .AND. (dnext .LT. dia)) then
							rate = dia / (dia - dnext)
							tout(kl) = tnow + rate * dt
						end if
						if (tout (kl) .GT.  now) then				! JMR: Whitespace!!!!!
							ddt = min(dt, (tout(kl) - tnow))
							wodot(kl) = (wo(kl) - wnext) / ddt
						else
							wodot(kl) = 0.0
						end if
						diam(kl) = dnext
						wo(kl) = wnext
					end if
					cycle lLoop
				end if

				!c If k of (k, l) still coming up to drying temperature, accumulate
				!c heat input and driving temperature difference, predict drying start

				if (tnow .LT. dryt) then
					factor = fmois(k) * dendry(k)
					conwet = condry(k) + 4.27e-04 * factor
					if (l .EQ. 0) then
						r = r0
						gi = fi + fid
					end if
					if (l .EQ. k) then
						r = r0
						gi = fi
					end if
					if ((l .NE. 0) .AND. (l .NE. k)) then
						r = r0 + 0.5 * flit(l) * dr
						gi = fi + flit(l) * fint(l)
					end if
					tf = tempf(gi, r, tpamb)
					if (tf .LE. (tpdry + 10.)) then
						cycle lLoop
					endif
					dia = diam(kl)
					ts = 0.5*(tpamb + tpdry)				! JMR: Whitespace!!!!!
					call heatx(u, d, dia, tf, ts, hf, hb, c, e)
					dtcum = max((tf - ts) * dt, 0.)
					tcum(kl) = tcum(kl) + dtcum
					qcum(kl) = qcum(kl) + hb * dtcum
					he = qcum(kl) / tcum(kl)
					dtef = tcum(kl) / tnext
					thd = (tpdry - tpamb) / dtef
					if (thd .GT. 0.9) then
						cycle lLoop
					endif
					biot = he * dia / conwet
					call DRYTIM(biot, thd, dryt)


! -- Pagebreak --
! Pg. 111:


					cpwet = cheat(k) + ch2o * fmois(k)
					fac = ((0.5 * dia) ** 2) / conwet
					fac = fac * cpwet * dendry(k)
					tdry(kl) = fac * dryt
					
					! JMR_TEMP_REPORTING:
					!if (kl .eq. 29) then ! Loc(7,1) = 29
					if (kl .eq. klVal) then
						print *, "tdry(kl) assignment 3: STEP()"
						print *, "tnow", tnow
						print *, "tdry(kl)", tdry(kl)
						print *, "fac", fac
					end if ! JMR_TEMP_REPORTING end
					
					if (tdry(kl) .LT. tnext) then
						ts = tpdry
						call heatx(u, d, dia, tf, ts, hf, hb, c, e)
						dqdt = hb * (tf - ts)
						delt = tnext - tdry(kl)
						qcum(kl) = dqdt * delt
						tcum(kl) = (tf - ts) * delt
						tbar = 0.5 * (tpdry + tpig(k))

						!c See if ignition to occur before time step complete

						if (tf .LE. (tpig(k) + 10.0)) then 
							cycle lLoop
						endif
						call TIGNIT(tpamb, tpdry, tpig(k), tf, &
									condry(k), cheat(k), fmois(k), dendry(k), &
									hb, dtlite)
						tign(kl) = 0.5 * (tdry(kl) + dtlite)
						
						! JMR_TEMP_REPORTING:
						!if (k .eq. 7) then ! Fuel type 7 is the first fuel showing issues.
						!if (kl .eq. 29) then! Loc(7,1) = 29
						if (kl .eq. klVal) then
							print *, "tign(kl) assignment 3: STEP()"
							print *, "tnow", tnow
							print *, "tign(kl)", tign(kl)
							!print *, "tpamb", tpamb
							!print *, "tpdry", tpdry
							!print *, "tpig(k)", tpig(k)
							print *, "tf", tf
							!print *, "condry(k)", condry(k)
							!print *, "cheat(k)", cheat(k)
							!print *, "fmois(k)", fmois(k)
							!print *, "dendry(k)", dendry(k)
							print *, "hb", hb
							print *, "dtlite", dtlite
							print *, "tdry(kl)", tdry(kl)
						end if ! JMR_TEMP_REPORTING end
						
						if (tnext .GT. tign(kl)) then
							ts = tchar(k)
							!qdot (kl, l) = hb * max ((tf - ts), 0.) ! JMR: Transcription: l for 1 substitution!!!!!
							qdot (kl, 1) = hb * max((tf - ts), 0.0)
						end if
					end if
				end if
			end do lLoop
		end do kLoop

		!c Update fractions ignited and burned out, to apply at next step start

		do k = 1, number
			flit (k) = 0.0
			fout (k) = 0.0
			do l = 0, k
				kl = Loc(k, l)
				flag = (tnext .GE. tign(kl))
				if (flag .AND. (tnext .LE. tout(kl))) then
					flit(k) = flit(k) + xmat(kl)
				end if
				if (tnext .GT. tout (kl)) then			! JMR_NOTE: Whitespace!!!!!
					fout(k) = fout(k) + xmat(kl)
				end if
			end do
		end do

	end subroutine STEP


! -- Pagebreak --
! Pg. 112:


	! This routine writes to file a summary of fuel consumption by size class, in one of two levels
	! of detail.
	!
	! History: Modernized original Burnup subroutine.
	! Note: While this reproduces the original output that output has some alignment issues.
	subroutine SUMMARY(outfil, number, maxno, maxkl, parts, nun, &
						tis, ak, wdry, fmois, sigma, tign, tout, xmat, wo, diam)
		implicit none

		! Arguments:
		character*12, intent(in) :: outfil			! Stores the name of input data files.
		! Note: outfil is defined right before the call to this routine and could be moved inside.
		integer, intent(in) :: number				! Actual number of fuel components
		integer, intent(in) :: maxno				! Max number of fuel components
		integer, intent(in) :: maxkl				! Max triangular matrix size.
		character*12, intent(in) :: parts(maxno)	! Fuel component names / labels
		integer, intent(in) :: nun					! Summary file unit identifier
		real, intent(in) :: tis						! Igniting surface fire residence time (s)		! JMR: This is a bit confusing.  The residence time (ti above) is output but tis is also the current time?
													! JMR_NOTE: No, this is definitely the 'current time', which when this is called will be the time the fire went out!!!!!
		real*4, intent(in) :: ak					! Area influence factor (ak / K_a parameter)
		real*4, intent(in) :: wdry(maxno)			! Ovendry mass loading, kg/sq m
		real*4, intent(in) :: fmois(maxno)			! Moisture content, fraction dry mass
		real*4, intent(in) :: sigma(maxno)			! Surface to volume ratio, 1 / m
		real*4, intent(in) :: tign(maxkl)			! Ignition time for the larger of each
													! fuel component pair
		real*4, intent(in) :: tout(maxkl)			! Burnout time of larger component of pairs
		real*4, intent(in) :: xmat(maxkl)			! Table of influence fractions between components
		real*4, intent(in) :: wo(maxkl)				! Current ovendry loading for the larger of
													! each component pair, kg/sq m
		real*4, intent(in) :: diam(maxkl)			! Current diameter of the larger of each
													! fuel component pair, m

		! Locals:
		character*3 :: stat		! File open status.
		integer :: ido, in		! User menu return value.
		logical :: v			! Should a full summary should be stored?
		character*12 :: nuname
		character*12 :: name
		integer :: m, n, mn		! Counter.
		real :: win, fmi, dim, rem, ts, tf
		real :: fr				! Value of xmat for a single element.
		real :: ti				! Value of tign for a single element.
		real :: t1				! Value of tout for a single element. (Was 'to', which is a keyword, in the original code.)
		real :: wd				! Value of wo for a single element.
		real :: di				! Diameter for a single element.
		integer :: openStat		! IO status.

		! Some of these format strings are so long that the 72 character limit renders them unreadable:
		! This format string is awkwardly long:
		character(len = *), parameter :: format20 = &
			"(/5x,a12,3x'Load , mois , diam =  '3e10.3/5x'  Companion  fraction   Ignition   burnout   load rem    diam')"
		character(len = *), parameter :: format30 = '(5x,a12,e9.3,e11.3,e10.3,e11.3,e9.3)'
		character(len = *), parameter :: format40 = '(5x,a12,e9.3,e11.3,e10.3,e11.3)'

		character*12, parameter :: noCmpStr = 'no companion'	! (Was a variable 'none', which is a
																!  keyword, in the original code.)
		
		nuname = outfil
		stat = 'NEW'

		! Ask for the summary type:
		do
			write(*, "(' Full summary = 1 , recap only = 0    ',$)")
			read(*, *) ido
			v = (ido .EQ. 1) ! JMR_NOTE: I think the use of 'v' makes the code less clear.
		!	if( ( ido .NE. 0) .AND. ( .NOT. v ) ) then
		!		cycle
		!	else
		!		exit
		!	end if
			! More compact version of the above code:
			if ((ido .eq. 0) .or. v) then
				exit
			end if
		end do

		! Open the summary file, assuming it does not exist yet:
		open(nun, file = nuname, status = stat, form = 'FORMATTED', iostat = openStat)

		! If there is an error opening the file prompt the user until resolved:
		do while (openStat .ne. 0)
			write(*, "('  Error on open:  'a12'  with status = 'a3)") nuname, stat
			write(*, *) ' 0 = abort now'
			write(*, *) ' 1 = overwrite existing file'!/
			write(*, "('  2 = enter new file name to create   ',$)")

			read(*,*) in
			! JMR_NOTE: Read error checking should be added here.
			if (in .EQ. 0) then
				stop ' Abort'
			else if (in .EQ. 1) then
				stat = 'OLD'
			else if (in .eq. 2) then
				write(*,*) 'Enter new file name [max 12 characters]'
				read(*, '(a12)') nuname
				stat = 'NEW'
			else
				cycle ! Ask again.
			end if
	
			! Will only get here if in == 1 or 2:
			open(nun, file = nuname, status = stat, form = 'FORMATTED', iostat = openStat)
		end do
		!goto 35
		! JMR_NOTE: This goto is a bug in the original code that will cause the first line, with the area
		! factor and end time, to be omitted from the summary file.  I am omitting it.

		! Output the summary:
		! JMR_Note: This is is very hard to write readably:
		write(nun, "(5x'Area factor =' , f6.2 , '  End time , s =' , e12.4 )") ak, tis

		if (.NOT. v) then ! ido must == 0, recap only
			write(nun, "(6x'Component   load @ 0   Ignition   burnout   load rem')")
		end if

		do m = 1, number
			name = parts(m)
			win = wdry(m)
			fmi = fmois(m)
			dim = 4. / sigma(m)
			if (v) then
				write(nun, format20) name, win, fmi, dim
			end if
	
			name = noCmpStr
			rem = 0.0
			ts = 1.0e+31
			tf = 0.0
			do n = 0, m
				mn = Loc(m, n)
				fr = xmat(mn)
				ti = tign(mn)
				ts = min(ti, ts)
				t1 = tout(mn)
				tf = max(t1, tf)
				wd = wo(mn)
				rem = rem + wd
				di = diam(mn)
				if (v) then
					write(nun, format30) name, fr, ti, t1, wd, di
				endif
				if (n .LT. m) then
					name = parts(n + 1)
				end if
			end do
			name = parts(m)
			if (.NOT. v) then
				write(nun, format40) name, win, ts, tf, rem
			end if


! -- Pagebreak --
! Pg. 113:


		end do

	end subroutine SUMMARY

! End of original source code.

	! This function coverts the indexes of the pairwise fuel interaction triangular matrix space
	! to indexes of the arrays used to represent it for several computed variables.
	!
	! Note: This will only return valid value valid (occupied) coordinates of the triangular matrix.
	! Error checking would require that the number of fuel classes be know. 
	!
	! History: This function was originally implemented as a statement function defined in seven places
	! in the original code.
	function Loc(k, l)
		 implicit none

		! Arguments:
		integer, intent(in) :: k	! Triangular matrix column (row) index, (1:number of fuel types)
		integer, intent(in) :: l	! Triangular matrix row (column) index, (0:k), = partner
									! This index starts at 0, which represent the "no companion"
									! pairs.

		! Locals:
		integer :: loc ! Return value: Index in a compact array representing the triangular matrix values.

		! Validity checking:
		if ((k .lt. 1) .or. (k .gt. maxno)) then
			print *, "Loc(): Invalid value of k ", k
		end if

		if ((k .lt. 0) .or. (k .gt. k)) then
			print *, "Loc(): Invalid value of l ", l
		end if

		loc = k * (k + 1) / 2 + l

		if ((loc .lt. 1) .or. (loc .gt. maxkl)) then
			print *, "Loc(): Invalid index returned ", loc
		end if

	end function Loc


	! I believe this function is the approximation of the complementary error function as
	! described in Albini 1995 and Albini & Reinhardt 1995.  It was obtained from Hastings (et
	! al.) 1955, but I cann't identify the equation in that reference.
	!
	! History: This function was originally implemented as a statement function in DRYTIM(), as f().
	function ErrorApprox(h, theta) result(approx)
		implicit none

		! Arguments:
		real*4, intent(in) :: h
		real*4, intent(in) :: theta	! Temperature rise required for the start of moisture loss.
	
		! Locals:
		real*4 :: approx ! Return value.
	
		! Constants:
		real*4, parameter :: a = 0.7478556
		real*4, parameter :: b = 0.4653628
		real*4, parameter :: c = 0.1282064

		approx = h * (b - h * (c - h)) - (1.0 - theta) / a

	end function ErrorApprox


	! AskForRealX() is a data entry utility function.  It prompts the user for a numeric value
	! (real) checking against the expected value range.  The value is returned.
	!
	! History: This function was implemented to reduce code repetition in GETDAT().
	! To reproduce the original behavior leading spaces should be omitted [in cName and paramDesc].
	function AskForReal1(cName, paramDesc, rangeLow, rangeHigh) result(input)
		implicit none

		! Arguments:
		character*12, intent(in) :: cName ! The fuel component name.
		character(len = *), intent(in) :: paramDesc ! A string describing the parameter.
		real, intent(in) :: rangeLow
		real, intent(in) :: rangeHigh

		! Constants:
		character(len = *), parameter :: fmtStart = "(1x,a12 '" ! Have to re-quote paramDesc.
		character(len = *), parameter :: fmtEnd = "',$)"

		! Locals:
		real :: input ! Return value.
		character(len = len(fmtStart) + len(paramDesc) + len(fmtEnd)) :: formatString
		integer :: readStat ! IO error status.
		logical :: startAgain

		formatString = fmtStart // paramDesc // fmtEnd

		do
			! Request the value:
			write(*, formatString) cName
		
			! Read the value:
			read(*, *, iostat = readStat) input
		
			if (readStat .ne. 0) then
				cycle ! If there was an error ask again.
			else
				! Check the value:
				call retry(input, rangeLow, rangeHigh, startAgain)
				if (startAgain) then
					cycle
				else
					exit
				end if
			end if
		end do

	end function AskForReal1


	! This version is for scalar paramaters and does not include a component name argument.
	function AskForReal2(paramDesc, rangeLow, rangeHigh) result(input)
		implicit none

		! Arguments:
		character(len = *), intent(in) :: paramDesc ! A string describing the parameter.
		real, intent(in) :: rangeLow
		real, intent(in) :: rangeHigh

		! Constants:
		character(len = *), parameter :: fmtStart = "(1x '" ! Have to re-quote paramDesc.
		character(len = *), parameter :: fmtEnd = "',$)"

		! Locals:
		real :: input ! Return value.
		character(len = len(fmtStart) + len(paramDesc) + len(fmtEnd)) :: formatString
		integer :: readStat ! IO error status.
		logical :: startAgain

		formatString = fmtStart // paramDesc // fmtEnd

		do
			! Request the value:
			write(*, trim(formatString))
		
			! Read the value:
			read(*, *, iostat = readStat) input
		
			if (readStat .ne. 0) then
				cycle ! If there was an error ask again.
			else
				! Check the value:
				call retry(input, rangeLow, rangeHigh, startAgain)
				if (startAgain) then
					cycle
				else
					exit
				end if
			end if
		end do

	end function AskForReal2


	! In this version the cName argument is optional.  This was intended to combine the behavior
	! of  allowing combines behavior of AskForReal1() and AskForReal2().  However, the compiler=
	! is giving errors when cName is omitted.
	function AskForReal(paramDesc, rangeLow, rangeHigh, cName) result(input)
		implicit none

		! Arguments:
		character(len = *), intent(in) :: paramDesc ! A string describing the parameter.
		real, intent(in) :: rangeLow ! The low end of the expected range.
		real, intent(in) :: rangeHigh ! The high end of the expected range.
		character*12, intent(in), optional :: cName ! Optional fuel component name.
	
		! Constants:
		character(len = *), parameter :: fmtStart = "(1x,a12 '" ! Have to re-quote paramDesc.
		character(len = *), parameter :: fmtEnd = "',$)"
 
		! Locals:
		real :: input ! Return value.
		character(len = len(fmtStart) + len(paramDesc) + len(fmtEnd)) :: formatString ! Longest case!
		integer :: readStat ! IO error status.
		logical :: startAgain
	
		! Compose the format string:
		if (present(cName)) then
			formatString = fmtStart // paramDesc // fmtEnd
		else
			formatString = "(1x, '" // paramDesc // fmtEnd
		end if
	
		do
			! Request the value:
			if (present(cName)) then
				write(*, formatString) cName
			else
				write(*, trim(formatString))
			end if
		
			! Read the value:
			read(*, *, iostat = readStat) input
		
			if (readStat .ne. 0) then
				cycle ! If there was an error ask again.
			else
				! Check the value:
				call retry(input, rangeLow, rangeHigh, startAgain)
				if (startAgain) then
					cycle
				else
					exit
				end if
			end if
		end do

	end function AskForReal

!end program BURNUP
end module BurnupMod

