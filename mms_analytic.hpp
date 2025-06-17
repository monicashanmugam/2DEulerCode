// mms_analytic.hpp
#ifndef MMS_ANALYTIC_HPP
#define MMS_ANALYTIC_HPP

#include <vector>
#include <cassert>
#include <cmath>
#include "variables.hpp"

// Declare the globals that were defined in main.cpp:
extern int imax;
extern int jmax;
extern std::vector<std::vector<Conserved>> U;
extern std::vector<std::vector<Primitive>> V;

//–– Prototypes of the four closed‐form “primitive” fields:
double rho_mms(int mmsCase, double L, double x, double y);
double uvel_mms(int mmsCase, double L, double x, double y);
double vvel_mms(int mmsCase, double L, double x, double y);
double press_mms(int mmsCase, double L, double x, double y);

//–– Prototypes of the four closed‐form “source” fields:
double rmassconv(int mmsCase, double L, double x, double y);
double xmtmconv(int mmsCase, double L, double x, double y);
double ymtmconv(int mmsCase, double L, double x, double y);
double energyconv(int mmsCase, double gamma, double L, double x, double y);

//–– Wrapper to initialize U/V to the exact MMS solution:
void initializeMMS(
  int mmsCase, 
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell,
  std::vector<std::vector<Conserved>>&   U,
  std::vector<std::vector<Primitive>>&   V
);

//–– Wrapper to compute the source‐term array S[i][j]:
void computeSourceTermsMMS(
  int mmsCase,
  double gamma,
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell,
  std::vector<std::vector<Conserved>>&   S
);

//---Wrapper to apply Dirichlet Boundary Conditions ---//

void applyMMSDirichletBCs(
  int mmsCase,
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell
);

void computeResidualMMS(
  int mmsCase,
  double gamma,
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell,
  const std::vector<std::vector<Primitive>>& V,
  std::vector<std::vector<Conserved>>& U,
  std::vector<std::vector<Conserved>>& R
);



void rungeKutta2Step(
  int mmsCase,
  double dt,
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell,
  const std::vector<std::vector<double>>& cellVolume,
  std::vector<std::vector<Conserved>>& U,
  std::vector<std::vector<Primitive>>& V
);



#endif // MMS_ANALYTIC_HPP
