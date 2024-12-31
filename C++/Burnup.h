/***************************************************************************************************
Burnup.h
Burnup Wildfire Fuel Consumption Model

Original Fortan code by: Frank A. Albini
Edited, modified, and ported to C++ by: Joshua M. Rady
Woodwell Climate Research Center
Started: 12/30/2024
Reference: Proj. 11 Exp. 22

This is an reimplementation of the Burnup wildfire fuel consumption model in C++.

...

***************************************************************************************************/
#ifndef BURNUP_H
#define BURNUP_H

//Core calculation functions:
//DUFBRN
//ARRAYS
//SORTER
//OVLAPS
//START
//FIRINT
//TIGNIT
double DRYTIM(const double enu, const double theta);
//HEATX
double TEMPF(const double q, const double r, const double tamb);
//STEP


//Utilities:
//Loc
double ErrorApprox(const double h, const double theta);
//AskForReal1
//AskForReal2
//AskForReal

#endif
