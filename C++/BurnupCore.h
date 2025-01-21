/***************************************************************************************************
BurnupCore.h
Burnup Wildfire Fuel Consumption Model

Original Fortan code by: Frank A. Albini
Edited, modified, and ported to C++ by: Joshua M. Rady
Woodwell Climate Research Center
Started: 12/30/2024
Reference: Proj. 11 Exp. 22

This is a reimplementation of the Burnup wildfire fuel consumption model in C++.
See Burnup.cpp and repository documentation for further inforation.

	There is no licence provided for the original code.  It is though to be open by provenance, but
that may not be correct.  The licence for this code is under consideration.
***************************************************************************************************/
#ifndef BURNUPCORE_H
#define BURNUPCORE_H

#include <string>

//Noninteractive Burnup entry point:
void Simulate(double& fi, const double ti, const double u, const double d, const double tpamb,
              const double ak, const double r0, const double dr, double& dt, const double wdf,
              const double dfm, const int ntimes, const int number,
              std::vector<std::string>& parts, std::vector<double>& wdry, std::vector<double>& ash,
              std::vector<double>& htval, std::vector<double>& fmois, std::vector<double>& dendry,
              std::vector<double>& sigma, std::vector<double>& cheat, std::vector<double>& condry,
              std::vector<double>& tpig, std::vector<double>& tchar, std::vector<double>& xmat,
              std::vector<double>& tign, std::vector<double>& tout, std::vector<double>& wo,
              std::vector<double>& diam, const bool outputHistory = false);

//Core calculation functions:
void DUFBRN(const double wdf, const double dfm, double& dfi, double& tdf);
void ARRAYS(std::vector<double>& wdry, std::vector<double>& ash, std::vector<double>& dendry,
            std::vector<double>& fmois, std::vector<double>& sigma, std::vector<double>& htval,
            std::vector<double>& cheat, std::vector<double>& condry, std::vector<double>& tpig,
            std::vector<double>& tchar, std::vector<double>& diam, std::vector<int>& key,
            std::vector<double>& work, const double ak, std::vector<std::vector<double>>& elam,
            std::vector<double>& alone, std::vector<double>& xmat, std::vector<double>& wo,
            std::vector<std::string>& parts, std::vector<std::string>& list,
            std::vector<double>& area, const int number = -1);
void SORTER(std::vector<double>& sigma, std::vector<double>& fmois, std::vector<double>& dryden,
            std::vector<int>& key, const int number = -1);
void OVLAPS(const std::vector<double> dryld, const std::vector<double> sigma,
            const std::vector<double> dryden, const double ak, const std::vector<double> fmois,
            std::vector<double>& beta, std::vector<std::vector<double>>& elam,
            std::vector<double>& alone, std::vector<double>& area, const int number = -1);
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
           const double r0, const double dr, int& ncalls, const int number = -1);
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
          const double tin, const std::vector<double> fint, const double fid, const int number = -1);

//Utilities:
int Loc(const int k, const int l);
int Length_kl(int numFuels);
double ErrorApprox(const double h, const double theta);
void ValidateOutputVector(std::vector<double>& output, const std:string outputName);

//File IO:
void SaveStateToFile(const int ts, const double time, const int number,
                     const std::vector<std::string> parts, const std::vector<double> wo,
                     const std::vector<double> diam, const double fi);

#endif
