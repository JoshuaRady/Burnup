# The Burnup Model, Old and New

This repository contains several variants of the Burnup wildfire fuel combustion model.  Developed in 1995 by Frank Albini Burnup is an influential computer model that simulates how much fuel will burn in a wildfire.  Burnup has subsequently been incorporated into other wildfire models such as FOFEM, FARSITE, and BehavePlus (Lutes 2013).  However, at the time this repository was created the source code for these models were not publicly accessible and the original model source code was hard to find.

This repository contains an annotated version of the original source code, which is in archaic fixed form Fortran, as well as a version that we have updated to modern Fortran.  Another version reimplemented as a Fortran module can be compiled as the original command line UI application or as a shared library with an additional interfaces that allows it to be run from R and has additional output options.

We have also ported the core calculations of Burnup to C++.  A R interface is included for this version as well.  Burnup has previously been ported to C++ but, as mentioned above, it is not publicly available.

All code versions give numerically identical results in testing.

This remains a work in progress with code verification ongoing.  Detailed documentation can be found in the code with Doxygen support in the C++ version.

## References

The orignal Burnup model was first documented in the report:

*Program BURNUP, a simulation model of the burning of large woody natural fuels.  
Frank A. Albini  
Missoula, MT: USDA Forest Service. Unpublished report. Research Grant INT-92754-GR. 1994.*

and the subsequent papers:

*Albini, F.A. and Reinhardt, E.D.  
Modeling ignition and burning rate of large woody natural fuels.  
International Journal of Wildland Fire 5(2): 81-91, 1995.*

*Albini, F.A., Brown, J.K., Reinhardt, E.D. and Ottmar, R.D.  
Calibration of a large fuel burnout model.  
International Journal of Wildland Fire 5(3): 173-192, 1995.*

This paper documents models that incorporate Burnup:

*D. C. Lutes.
Predicted fuel consumption in the Burnup model: sensitivity to four user inputs.
Res. Note RMRS-RN-51WWW. Fort Collins, CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. 11 p. https://doi.org/10.2737/RMRS-RN-51*
