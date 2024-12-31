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
void FIRINT(const std::vector<double> wodot, const std::vector<double> ash,
            const std::vector<double> htval, const int number, const std::vector<double> area,
            std::vector<double>& fint, double& fi);
double TIGNIT(const double tpam, const double tpdr, const double tpig, const double tpfi,
              const double cond, const double chtd, const double fmof, const double dend,
              const double hbar);
double DRYTIM(const double enu, const double theta);
void HEATX(const double u, const double d, const double dia, const double tf, const double ts,
           double& hfm, double& hbar, const double cond, double& en);
double TEMPF(const double q, const double r, const double tamb);
//STEP


//Utilities:
int Loc(const int k, const int l);
double ErrorApprox(const double h, const double theta);
//AskForReal1
//AskForReal2
//AskForReal

#endif
