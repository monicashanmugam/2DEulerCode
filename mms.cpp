#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <physicalContants.hpp>

// Supersonic Case : Table A.1 constants

namespace Supersonic {
  constexpr double rho0   = 1.0;
  constexpr double rhox   =  0.15;  // linear‐gradient coefficient for ρ
  constexpr double rhoy   = -0.10;  // linear‐gradient coefficient for ρ
  constexpr double arhox  =  1.0;   // amplitude for sin(π x/L)
  constexpr double arhoy  =  0.5;   // amplitude for sin(π y/L)

  constexpr double u0     = 800.0;
  constexpr double ux     =  50.0;  // linear‐gradient for u
  constexpr double uy     = -30.0;
  constexpr double aux    =   1.5;  // amplitude for sin(π x/L)
  constexpr double auy    =   0.6;  // amplitude for sin(π y/L)

  constexpr double v0     = 800.0;
  constexpr double vx     = -75.0;  // linear‐gradient for v
  constexpr double vy     =  40.0;
  constexpr double avx    =   0.5;  // amplitude for cos(π x/L)
  constexpr double avy    = 2.0/3.0;// amplitude for sin(π y/L)

  constexpr double p0     = 1.0e5;
  constexpr double px     = 0.2e5; // linear‐gradient for p
  constexpr double py     = 0.5e5;
  constexpr double apx    =   2.0; // amplitude for cos(π x/L)
  constexpr double apy    =   1.0; // amplitude for sin(π y/L)
}

// Subsonic Case : Table A.2 constants

namespace Subsonic {
  constexpr double rho0   = 1.0;
  constexpr double rhox   =  0.15;
  constexpr double rhoy   = -0.10;
  constexpr double arhox  =  1.0;
  constexpr double arhoy  =  0.5;

  constexpr double u0     =  70.0;
  constexpr double ux     =   5.0;
  constexpr double uy     =  -7.0;
  constexpr double aux    =   1.5;
  constexpr double auy    =   0.6;

  constexpr double v0     =  90.0;
  constexpr double vx     = -15.0;
  constexpr double vy     =   8.5;
  constexpr double avx    =   0.5;
  constexpr double avy    = 2.0/3.0;

  constexpr double p0     = 1.0e5;
  constexpr double px     = 0.2e5;
  constexpr double py     = 0.5e5;
  constexpr double apx    =   2.0;
  constexpr double apy    =   1.0;
}

struct Conserved {
  double rho;    // density
  double rhou;   // ρ·u
  double rhov;   // ρ·v
  double E;      // ρ·e_t  (total energy)
};

// Helper: init supersonic/subsonic fills entire array of the conserved variable U at a single point x,y

static void initMMS_supersonic(double x, double y,
                               double &rho, double &u, double &v,
                               double &p,   double &T, double &et)
{
  using namespace Supersonic;
  const double pi = M_PI;

  // 1) density = linear + sinusoidal
  rho = rho0
      + (rhox * (x / L))
      + (rhoy * (y / L))
      + (arhox * std::sin(pi * x / L))
      + (arhoy * std::sin(pi * y / L));

  // 2) u‐velocity = linear + sinusoidal
  u   = u0
      + (ux * (x / L))
      + (uy * (y / L))
      + (aux * std::sin(pi * x / L))
      + (auy * std::sin(pi * y / L));

  // 3) v‐velocity = linear + sinusoidal
  v   = v0
      + (vx * (x / L))
      + (vy * (y / L))
      + (avx * std::cos(pi * x / L))
      + (avy * std::sin(pi * y / L));

  // 4) pressure = linear + sinusoidal
  p   = p0
      + (px * (x / L))
      + (py * (y / L))
      + (apx * std::cos(pi * x / L))
      + (apy * std::sin(pi * y / L));

  // 5) temperature via ideal‐gas law
  T   = p / (rho * Rgas);

  // 6) total energy per unit mass e_t = internal + kinetic
  double e_internal = p / ((gamma_gas - 1.0) * rho);
  et   = e_internal + 0.5*(u*u + v*v);
}

static void initMMS_subsonic(double x, double y,
                             double &rho, double &u, double &v,
                             double &p,   double &T, double &et)
{
  using namespace Subsonic;
  const double pi = M_PI;

  // 1) density
  rho = rho0
      + (rhox * (x / L))
      + (rhoy * (y / L))
      + (arhox * std::sin(pi * x / L))
      + (arhoy * std::sin(pi * y / L));

  // 2) u‐velocity
  u   = u0
      + (ux * (x / L))
      + (uy * (y / L))
      + (aux * std::sin(pi * x / L))
      + (auy * std::sin(pi * y / L));

  // 3) v‐velocity
  v   = v0
      + (vx * (x / L))
      + (vy * (y / L))
      + (avx * std::cos(pi * x / L))
      + (avy * std::sin(pi * y / L));

  // 4) pressure
  p   = p0
      + (px * (x / L))
      + (py * (y / L))
      + (apx * std::cos(pi * x / L))
      + (apy * std::sin(pi * y / L));

  // 5) temperature
  T   = p / (rho * Rgas);

  // 6) total energy
  double e_internal = p / ((gamma_gas - 1.0) * rho);
  et   = e_internal + 0.5*(u*u + v*v);
}

// The wrapper : InitializeMMS() Fills the entire 2D space U[i][j] with exact MMS conserved state
// "mmsCase" lets you choose between 1) supersonic 2) subsonic

static void InitializeMMS(
    const std::vector<std::vector<double>>& xcc,
    const std::vector<std::vector<double>>& ycc,
    std::vector<std::vector<Conserved>>&   U,
    int                                     mmsCase)
{
    assert(mmsCase == 1 || mmsCase == 2);

    int Ni = static_cast<int>(xcc.size());      // = imax + 2*ghost
    int Nj = static_cast<int>(xcc[0].size());   // = jmax + 2*ghost

    // Resize U if needed
    U.assign(Ni, std::vector<Conserved>(Nj));

    // Loop over every (i,j), including ghost layers if you like.
    for (int i = 0; i < Ni; ++i) {
      for (int j = 0; j < Nj; ++j) {
        double x = xcc[i][j];
        double y = ycc[i][j];

        double rho, u, v, p, T, et;
        if (mmsCase == 1) {
          initMMS_supersonic(x, y, rho, u, v, p, T, et);
        } else {
          initMMS_subsonic(x, y, rho, u, v, p, T, et);
        }

        // Store into 2D conserved structure (2D Euler, so no w‐component)
        U[i][j].rho  = rho;
        U[i][j].rhou = rho * u;
        U[i][j].rhov = rho * v;
        U[i][j].E    = rho * et;
      }
    }
}
