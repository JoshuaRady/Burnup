!---------------------------------------------------------------------------------------------------
! BurnupInteractive.f90
! Burnup Wildfire Fuel Consumption Model
!
! Programmed by : Joshua M. Rady
! Woodwell Climate Research Center
! Started: 12/19/2023
!
! This is a wrapper for compiling the modernized Burnup module as an interactive program with the
! menu driven command line interface of the original Burnup implementation (Albini 1994).
!
!---------------------------------------------------------------------------------------------------

program BurnupInteractive

	use BurnupMod, only : InteractiveUI

	implicit none

	call InteractiveUI() ! Call the interactive UI entry point.

end program BurnupInteractive
