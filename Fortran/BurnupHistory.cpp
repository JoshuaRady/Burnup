/***************************************************************************************************
BurnupHistory.cpp
Programmed by: Joshua M. Rady
Woodwell Climate Research Center
Started: 7/18/2025
Reference: Proj. 11 Exp. 25

	This provides an object to store Burnup model state during the time evolution of a simulation.

Licence?????
***************************************************************************************************/

#include "BurnupHistory.h"

//We need a persistant object to short data to.  This is kept private in this file, only being\
//accessed via the provided functions:
BurnupHistory BUHistStore;

void SaveStateToHistory(const int ts, const double time, const int number,
                        const std::vector<std::string> parts, const std::vector<double> wo,
                        const std::vector<double> diam, const double fi)
{
	
}


BurnupHistory GetHistory()
{
	//Add checking that the history is complete?
	return BUHistStore;
}