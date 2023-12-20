!---------------------------------------------------------------------------------------------------
! BurnupInteractive.f90
!
! Programmed by : Joshua M. Rady
! Started: 12/19/2023
!
! This is a wrapper for compiling the Burnup module as an interactive program, similar to the
! original implementation.
!
!---------------------------------------------------------------------------------------------------

program Wrapper

	use BurnupMod, only : BurnupMain

	implicit none

	call BurnupMain()

end program Wrapper
