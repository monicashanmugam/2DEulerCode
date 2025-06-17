// mms_analytic.cpp
#include "mms_analytic.hpp"
#include "variables.hpp"
#include "physicalConstants.hpp"

// Declare the globals that were defined in main.cpp:
extern int imax;
extern int jmax;
extern std::vector<std::vector<Conserved>> U;
extern std::vector<std::vector<Primitive>> V;

//----------------------------------------------------------------------
// (A) Define supersonic (Table A.1) and subsonic (Table A.2) constants
//----------------------------------------------------------------------
namespace {
  // Supersonic constants (Table A.1)
  struct MmsParams {
    double rho0,   rho_x,   rho_y,   a_rho_x, a_rho_y;
    double u0,     u_x,     u_y,     a_u_x,   a_u_y;
    double v0,     v_x,     v_y,     a_v_x,   a_v_y;
    double p0,     p_x,     p_y,     a_p_x,   a_p_y;
  };

  static constexpr MmsParams mmsSup = {
    1.0,   0.15,  -0.10,  1.0,    0.50,    // rho…
    800.0, 50.0, -30.0,  1.5,    0.60,    // u…
    800.0, -75.0, 40.0,  0.5,    (2.0/3.0), // v…
    100000.0, 20000.0, 50000.0, 2.0, 1.0   // p…
  };

  // Subsonic constants (Table A.2)
  static constexpr MmsParams mmsSub = {
    1.0,    0.15,  -0.10,  1.0,   0.50,   // rho…
     70.0,   5.0,   -7.0,  1.5,  0.60,    // u…
     90.0,  -15.0,   8.5,  0.5,  (2.0/3.0), // v…
   100000.0, 20000.0, 50000.0, 2.0, 1.0   // p…
  };

  // A shared π constant (you may also use M_PI if <cmath> provides it):
  static constexpr double PI = std::acos(-1.0);
}

//----------------------------------------------------------------------
// (B) Exact “primitive” fields (density, u, v, pressure)
//----------------------------------------------------------------------
double rho_mms(int mmsCase, double L, double x, double y) {
  assert(mmsCase==1 || mmsCase==2);
  auto const &C = (mmsCase==1 ? mmsSup : mmsSub);

  double ξ   = x / L;
  double η   = y / L;
  double lin = C.rho0 + C.rho_x*ξ + C.rho_y*η;
  double sin = C.a_rho_x*std::sin(PI*ξ)
             + C.a_rho_y*std::sin(PI*η);
  return lin + sin;
}

double uvel_mms(int mmsCase, double L, double x, double y) {
  assert(mmsCase==1 || mmsCase==2);
  auto const &C = (mmsCase==1 ? mmsSup : mmsSub);

  double ξ   = x / L;
  double η   = y / L;
  double lin = C.u0 + C.u_x*ξ + C.u_y*η;
  double sin = C.a_u_x*std::sin(PI*ξ)
             + C.a_u_y*std::sin(PI*η);
  return lin + sin;
}

double vvel_mms(int mmsCase, double L, double x, double y) {
  assert(mmsCase==1 || mmsCase==2);
  auto const &C = (mmsCase==1 ? mmsSup : mmsSub);

  double ξ   = x / L;
  double η   = y / L;
  double lin = C.v0 + C.v_x*ξ + C.v_y*η;
  double sin = C.a_v_x*std::sin(PI*ξ)
             + C.a_v_y*std::sin(PI*η);
  return lin + sin;
}

double press_mms(int mmsCase, double L, double x, double y) {
  assert(mmsCase==1 || mmsCase==2);
  auto const &C = (mmsCase==1 ? mmsSup : mmsSub);

  double ξ   = x / L;
  double η   = y / L;
  double lin = C.p0 + C.p_x*ξ + C.p_y*η;
  double sin = C.a_p_x*std::sin(PI*ξ)
             + C.a_p_y*std::sin(PI*η);
  return lin + sin;
}

//----------------------------------------------------------------------
// (C) Exact “source” fields (mass‐equation, x‐momentum, y‐momentum, energy)
//   Each formula is copied verbatim from your Fortran → C++ conversion,
//   but replacing every “rho0, rhox, …” with C.<name> based on mmsCase.
//----------------------------------------------------------------------

double rmassconv(int mmsCase, double L, double x, double y) {
  assert(mmsCase==1 || mmsCase==2);
  auto const &C = (mmsCase==1 ? mmsSup : mmsSub);

  // term1
  double t1 = (3.0*PI*C.u_x * std::cos(3.0*PI*x/(2.0*L))
             * (C.rho0 + C.rho_y*std::cos(PI*y/(2.0*L))
                       + C.rho_x*std::sin(PI*x/L)))
            / (2.0*L);

  // term2
  double t2 = (2.0*PI*C.v_y * std::cos(2.0*PI*y/(3.0*L))
             * (C.rho0 + C.rho_y*std::cos(PI*y/(2.0*L))
                       + C.rho_x*std::sin(PI*x/L)))
            / (3.0*L);

  // term3
  double t3 = (PI*C.rho_x * std::cos(PI*x/L)
             * (C.u0 + C.u_y*std::cos(3.0*PI*y/(5.0*L))
                     + C.u_x*std::sin(3.0*PI*x/(2.0*L))))
            / L;

  // term4
  double t4 = (PI*C.rho_y * std::sin(PI*y/(2.0*L))
             * (C.v0 + C.v_x*std::cos(PI*x/(2.0*L))
                       + C.v_y*std::sin(2.0*PI*y/(3.0*L))))
            / (2.0*L);

  return t1 + t2 + t3 - t4;
}

double xmtmconv(int mmsCase, double L, double x, double y) {
  assert(mmsCase==1 || mmsCase==2);
  auto const &C = (mmsCase==1 ? mmsSup : mmsSub);

  using std::sin;
  using std::cos;
  using std::pow;

  double t1 = (3.0*PI*C.u_x * cos(3.0*PI*x/(2.0*L))
             * (C.rho0 + C.rho_y*cos(PI*y/(2.0*L))
                       + C.rho_x*sin(PI*x/L))
             * (C.u0 + C.u_y*cos(3.0*PI*y/(5.0*L))
                       + C.u_x*sin(3.0*PI*x/(2.0*L))))
            / L;

  double t2 = (2.0*PI*C.v_y * cos(2.0*PI*y/(3.0*L))
             * (C.rho0 + C.rho_y*cos(PI*y/(2.0*L))
                       + C.rho_x*sin(PI*x/L))
             * (C.u0 + C.u_y*cos(3.0*PI*y/(5.0*L))
                       + C.u_x*sin(3.0*PI*x/(2.0*L))))
            / (3.0*L);

  double t3 = (PI*C.rho_x * cos(PI*x/L)
             * pow(C.u0 + C.u_y*cos(3.0*PI*y/(5.0*L))
                     + C.u_x*sin(3.0*PI*x/(2.0*L)), 2.0))
            / L;

  double t4 = -(2.0*PI*C.p_x * sin(2.0*PI*x/L)) / L;

  double t5 = -(PI*C.rho_y
              * (C.u0 + C.u_y*cos(3.0*PI*y/(5.0*L))
                        + C.u_x*sin(3.0*PI*x/(2.0*L)))
              * sin(PI*y/(2.0*L))
              * (C.v0 + C.v_x*cos(PI*x/(2.0*L))
                        + C.v_y*sin(2.0*PI*y/(3.0*L))))
            / (2.0*L);

  return t1 + t2 + t3 + t4 + t5;
}

double ymtmconv(int mmsCase, double L, double x, double y) {
  assert(mmsCase==1 || mmsCase==2);
  auto const &C = (mmsCase==1 ? mmsSup : mmsSub);

  using std::sin;
  using std::cos;
  using std::pow;

  double t1 = (PI * C.p_y * cos(PI*y/L)) / L;

  double t2 = -(PI * C.v_x * sin(PI*x/(2.0*L))
             * (C.rho0 + C.rho_y * cos(PI*y/(2.0*L))
                       + C.rho_x * sin(PI*x/L))
             * (C.u0 + C.u_y * cos(3.0*PI*y/(5.0*L))
                       + C.u_x * sin(3.0*PI*x/(2.0*L))))
            / (2.0*L);

  double t3 = (3.0*PI * C.u_x * cos(3.0*PI*x/(2.0*L))
             * (C.rho0 + C.rho_y * cos(PI*y/(2.0*L))
                       + C.rho_x * sin(PI*x/L))
             * (C.v0 + C.v_x * cos(PI*x/(2.0*L))
                       + C.v_y * sin(2.0*PI*y/(3.0*L))))
            / (2.0*L);

  double t4 = (4.0*PI * C.v_y * cos(2.0*PI*y/(3.0*L))
             * (C.rho0 + C.rho_y * cos(PI*y/(2.0*L))
                       + C.rho_x * sin(PI*x/L))
             * (C.v0 + C.v_x * cos(PI*x/(2.0*L))
                       + C.v_y * sin(2.0*PI*y/(3.0*L))))
            / (3.0*L);

  double t5 = (PI * C.rho_x * cos(PI*x/L)
             * (C.u0 + C.u_y * cos(3.0*PI*y/(5.0*L))
                       + C.u_x * sin(3.0*PI*x/(2.0*L)))
             * (C.v0 + C.v_x * cos(PI*x/(2.0*L))
                       + C.v_y * sin(2.0*PI*y/(3.0*L))))
            / L;

  double t6 = -(PI * C.rho_y * sin(PI*y/(2.0*L))
             * pow(C.v0 + C.v_x * cos(PI*x/(2.0*L))
                   + C.v_y * sin(2.0*PI*y/(3.0*L)), 2.0))
            / (2.0*L);

  return t1 + t2 + t3 + t4 + t5 + t6;
}

double energyconv(int mmsCase, double gamma, double L, double x, double y) {
  assert(mmsCase==1 || mmsCase==2);
  auto const &C = (mmsCase==1 ? mmsSup : mmsSub);

  using std::sin;
  using std::cos;
  using std::pow;

  // term1 = (u…) * [ … ]  –  (π·rhox·cos…)/(…)
  double term1 = ( C.u0
                  + C.u_y * cos((3.0*PI*y)/(5.0*L))
                  + C.u_x * sin((3.0*PI*x)/(2.0*L))
                 ) * (
                   (-2.0*PI*C.p_x * sin((2.0*PI*x)/L)) / L
                   + ( C.rho0
                     + C.rho_y * cos((PI*y)/(2.0*L))
                     + C.rho_x * sin((PI*x)/L)
                     ) * (
                       (-2.0*PI*C.p_x * sin((2.0*PI*x)/L))
                         / ((-1.0 + gamma)*L
                            * (C.rho0
                              + C.rho_y * cos((PI*y)/(2.0*L))
                              + C.rho_x * sin((PI*x)/L)
                              )
                            )
                       + (
                           (3.0*PI*C.u_x * cos((3.0*PI*x)/(2.0*L))
                            * ( C.u0
                              + C.u_y * cos((3.0*PI*y)/(5.0*L))
                              + C.u_x * sin((3.0*PI*x)/(2.0*L))
                              )
                           )/L
                           - (PI*C.v_x * sin((PI*x)/(2.0*L))
                              * ( C.v0
                                + C.v_x * cos((PI*x)/(2.0*L))
                                + C.v_y * sin((2.0*PI*y)/(3.0*L))
                                )
                              )/L
                         )/2.0
                     )
                 )
                 - (PI*C.rho_x * cos((PI*x)/L)
                    * ( C.p0
                      + C.p_x * cos((2.0*PI*x)/L)
                      + C.p_y * sin((PI*y)/L)
                      )
                   ) / ((-1.0 + gamma)*L
                        * pow(C.rho0
                              + C.rho_y * cos((PI*y)/(2.0*L))
                              + C.rho_x * sin((PI*x)/L)
                              , 2.0)
                         );

  // term2 = (π·rhox·cos…)*[ … ]/L
  double term2 = ( PI*C.rho_x * cos((PI*x)/L) * (
                   ( pow(C.u0
                          + C.u_y * cos((3.0*PI*y)/(5.0*L))
                          + C.u_x * sin((3.0*PI*x)/(2.0*L))
                        , 2.0)
                     + pow(C.v0
                           + C.v_x * cos((PI*x)/(2.0*L))
                           + C.v_y * sin((2.0*PI*y)/(3.0*L))
                         , 2.0)
                   )/2.0
                   + ( C.p0
                     + C.p_x * cos((2.0*PI*x)/L)
                     + C.p_y * sin((PI*y)/L)
                     ) / ((-1.0 + gamma)
                          * (C.rho0
                            + C.rho_y * cos((PI*y)/(2.0*L))
                            + C.rho_x * sin((PI*x)/L)
                            )
                          )
                 )) / L;

  // term3 = (3·π·uvelx·cos…)*[ … ]/(2L)
  double term3 = ( 3.0*PI*C.u_x * cos((3.0*PI*x)/(2.0*L)) * (
                    (C.p0
                     + C.p_x * cos((2.0*PI*x)/L)
                     + C.p_y * sin((PI*y)/L)
                    )
                    + ( C.rho0
                      + C.rho_y * cos((PI*y)/(2.0*L))
                      + C.rho_x * sin((PI*x)/L)
                      ) * (
                        ( pow(C.u0
                               + C.u_y * cos((3.0*PI*y)/(5.0*L))
                               + C.u_x * sin((3.0*PI*x)/(2.0*L))
                             , 2.0)
                          + pow(C.v0
                                + C.v_x * cos((PI*x)/(2.0*L))
                                + C.v_y * sin((2.0*PI*y)/(3.0*L))
                              , 2.0)
                        )/2.0
                        + (C.p0
                           + C.p_x * cos((2.0*PI*x)/L)
                           + C.p_y * sin((PI*y)/L)
                          )/(( -1.0 + gamma)
                              * (C.rho0
                                + C.rho_y * cos((PI*y)/(2.0*L))
                                + C.rho_x * sin((PI*x)/L)
                                )
                              )
                      )
                  )) / (2.0*L);

  // term4 = (2·π·vvely·cos…)*[ … ]/(3L)
  double term4 = ( 2.0*PI*C.v_y * cos((2.0*PI*y)/(3.0*L)) * (
                    (C.p0
                     + C.p_x * cos((2.0*PI*x)/L)
                     + C.p_y * sin((PI*y)/L)
                    )
                    + ( C.rho0
                      + C.rho_y * cos((PI*y)/(2.0*L))
                      + C.rho_x * sin((PI*x)/L)
                      ) * (
                        ( pow(C.u0
                               + C.u_y * cos((3.0*PI*y)/(5.0*L))
                               + C.u_x * sin((3.0*PI*x)/(2.0*L))
                             , 2.0)
                          + pow(C.v0
                                + C.v_x * cos((PI*x)/(2.0*L))
                                + C.v_y * sin((2.0*PI*y)/(3.0*L))
                              , 2.0)
                        )/2.0
                        + (C.p0
                           + C.p_x * cos((2.0*PI*x)/L)
                           + C.p_y * sin((PI*y)/L)
                          )/(( -1.0 + gamma)
                              * (C.rho0
                                + C.rho_y * cos((PI*y)/(2.0*L))
                                + C.rho_x * sin((PI*x)/L)
                                )
                              )
                      )
                  )) / (3.0*L);

  // term5 = (v…) * [dG/dy part] / (2L)
  double term5 = ( C.v0
                  + C.v_x * cos((PI*x)/(2.0*L))
                  + C.v_y * sin((2.0*PI*y)/(3.0*L))
                 ) * (
                   (PI*C.p_y * cos((PI*y)/L)) / L
                   - (PI*C.rho_y * sin((PI*y)/(2.0*L))
                      * (
                          ( pow(C.u0
                                 + C.u_y * cos((3.0*PI*y)/(5.0*L))
                                 + C.u_x * sin((3.0*PI*x)/(2.0*L))
                               , 2.0)
                            + pow(C.v0
                                  + C.v_x * cos((PI*x)/(2.0*L))
                                  + C.v_y * sin((2.0*PI*y)/(3.0*L))
                                 , 2.0)
                          )/2.0
                          + (C.p0
                             + C.p_x * cos((2.0*PI*x)/L)
                             + C.p_y * sin((PI*y)/L)
                            )/(( -1.0 + gamma)
                                * (C.rho0
                                  + C.rho_y * cos((PI*y)/(2.0*L))
                                  + C.rho_x * sin((PI*x)/L)
                                  )
                                )
                        )
                     ) / (2.0*L)
                 );

  // term6 = (rho…) * [combination of derivatives] 
  double term6 = ( C.rho0
                  + C.rho_y * cos((PI*y)/(2.0*L))
                  + C.rho_x * sin((PI*x)/L)
                 ) * (
                   (PI*C.p_y * cos((PI*y)/L))
                     / ((-1.0 + gamma)*L
                        * (C.rho0
                          + C.rho_y * cos((PI*y)/(2.0*L))
                          + C.rho_x * sin((PI*x)/L)
                          )
                        )
                   + (
                       (-6.0*PI*C.u_y
                         * (C.u0
                            + C.u_y * cos((3.0*PI*y)/(5.0*L))
                            + C.u_x * sin((3.0*PI*x)/(2.0*L))
                           )
                         * sin((3.0*PI*y)/(5.0*L))
                       ) / (5.0*L)
                       + (4.0*PI*C.v_y
                          * cos((2.0*PI*y)/(3.0*L))
                          * (C.v0
                             + C.v_x * cos((PI*x)/(2.0*L))
                             + C.v_y * sin((2.0*PI*y)/(3.0*L))
                            )
                         ) / (3.0*L)
                     ) / 2.0
                   + (PI*C.rho_y * sin((PI*y)/(2.0*L))
                      * (C.p0
                         + C.p_x * cos((2.0*PI*x)/L)
                         + C.p_y * sin((PI*y)/L)
                        )
                     ) / (2.0 * (-1.0 + gamma) * L
                            * pow(C.rho0
                                  + C.rho_y * cos((PI*y)/(2.0*L))
                                  + C.rho_x * sin((PI*x)/L)
                                  , 2.0)
                           )
                 );

  return term1 + term2 + term3 + term4 + term5 + term6;
}

//----------------------------------------------------------------------
// (D) initializeMMS(…)
//----------------------------------------------------------------------
void initializeMMS(
  int mmsCase, 
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell,
  std::vector<std::vector<Conserved>>&   U,
  std::vector<std::vector<Primitive>>&   V
) {
  assert(mmsCase==1 || mmsCase==2);

  // Resize (including ghosts):
  int Ni = imax + 2*ghost;
  int Nj = jmax + 2*ghost;
  U.assign(Ni, std::vector<Conserved>(Nj));
  V.assign(Ni, std::vector<Primitive>(Nj));

  // Fill interior cells only:
  for(int i = ghost; i < ghost + imax; ++i) {
    for(int j = ghost; j < ghost + jmax; ++j) {
      double xc = x_cell[i][j];
      double yc = y_cell[i][j];

      double rho_ex = rho_mms   (mmsCase, L, xc, yc);
      double  u_ex = uvel_mms (mmsCase, L, xc, yc);
      double  v_ex = vvel_mms (mmsCase, L, xc, yc);
      double  p_ex = press_mms(mmsCase, L, xc, yc);

      V[i][j] = Primitive{rho_ex, u_ex, v_ex, p_ex};
      U[i][j] = PrimitiveToConserved(V[i][j]);
    }
  }
  // Ghost cells (i<ghost or i>=ghost+imax, or j<ghost or j>=ghost+jmax)
  // remain whatever your solver normally uses. (Often left at zero, or copied
  // from nearest interior, depending on boundary‐condition strategy.)
}

//----------------------------------------------------------------------
// (E) computeSourceTermsMMS(…)
//----------------------------------------------------------------------
void computeSourceTermsMMS(
  int mmsCase,
  double gamma,
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell,
  std::vector<std::vector<Conserved>>&   S
) {
  assert(mmsCase==1 || mmsCase==2);

  // Allocate (including ghosts) and zero everything:
  int Ni = imax + 2*ghost;
  int Nj = jmax + 2*ghost;
  S.assign(Ni, std::vector<Conserved>(Nj, Conserved{0,0,0,0}));

  // Fill interior cells only:
  for(int i = ghost; i < ghost + imax; ++i) {
    for(int j = ghost; j < ghost + jmax; ++j) {
      double xc = x_cell[i][j];
      double yc = y_cell[i][j];

      S[i][j].rho  = rmassconv   (mmsCase, L, xc, yc);
      S[i][j].rhou = xmtmconv    (mmsCase, L, xc, yc);
      S[i][j].rhov = ymtmconv    (mmsCase, L, xc, yc);
      S[i][j].E    = energyconv  (mmsCase, gamma, L, xc, yc);
    }
  }
  // Ghost rows of S[] remain zero—no source outside the interior.
}

//-----------------------------------------//
//Setting the  Dirichlet Boundary conditions
//-----------------------------------------//

void applyMMSDirichletBCs(
  int mmsCase,
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell
) {
  // Apply BCs to bottom and top ghost cells (j = 0,1 and jmax+2, jmax+3)
  for (int i = 0; i <= imax + 3; ++i) {
    for (int j : {0, 1, jmax + 2, jmax + 3}) {
      double x = x_cell[i][j];
      double y = y_cell[i][j];
      double rho = rho_mms(mmsCase, L, x, y);
      double u   = uvel_mms(mmsCase, L, x, y);
      double v   = vvel_mms(mmsCase, L, x, y);
      double p   = press_mms(mmsCase, L, x, y);
      V[i][j] = {rho, u, v, p};
      U[i][j] = primitive_to_conserved(V[i][j]);
    }
  }

  // Apply BCs to left and right ghost cells (i = 0,1 and imax+2, imax+3)
  for (int j = 0; j <= jmax + 3; ++j) {
    for (int i : {0, 1, imax + 2, imax + 3}) {
      double x = x_cell[i][j];
      double y = y_cell[i][j];
      double rho = rho_mms(mmsCase, L, x, y);
      double u   = uvel_mms(mmsCase, L, x, y);
      double v   = vvel_mms(mmsCase, L, x, y);
      double p   = press_mms(mmsCase, L, x, y);
      V[i][j] = {rho, u, v, p};
      U[i][j] = primitive_to_conserved(V[i][j]);
    }
  }
}

void rungeKutta2Step(
  int mmsCase,
  double dt,
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell,
  const std::vector<std::vector<double>>& cellVolume,
  std::vector<std::vector<Conserved>>& U,
  std::vector<std::vector<Primitive>>& V
) {
  // 1. Stage 1 residual
  std::vector<std::vector<Conserved>> R(imax+4, std::vector<Conserved>(jmax+4));
  computeResidualMMS(mmsCase, gamma_, L, x_cell, y_cell, V, U, R);

  std::vector<std::vector<Conserved>> U_stage(imax+4, std::vector<Conserved>(jmax+4));
  for (int i = 2; i <= imax+1; ++i) {
    for (int j = 2; j <= jmax+1; ++j) {
      double vol = cellVolume[i][j];
      U_stage[i][j] = U[i][j] - (0.5 * dt / vol) * R[i][j];
    }
  }
  // Enforce BC on U_stage
  U = U_stage;
  for (int i = 2; i <= imax+1; ++i)
    for (int j = 2; j <= jmax+1; ++j)
      V[i][j] = conserved_to_primitive(U[i][j]);

  applyMMSDirichletBCs(mmsCase, L, x_cell, y_cell);

  // 2. Stage 2 residual
  std::vector<std::vector<Conserved>> R2(imax+4, std::vector<Conserved>(jmax+4));
 computeResidualMMS(mmsCase, gamma_, L, x_cell, y_cell, V, U, R);

  for (int i = 2; i <= imax+1; ++i) {
    for (int j = 2; j <= jmax+1; ++j) {
      double vol = cellVolume[i][j];
      U[i][j] = U[i][j] - (1.0 * dt / vol) * R2[i][j];
    }
  }

  // Final BC and primitive update
  for (int i = 2; i <= imax+1; ++i)
    for (int j = 2; j <= jmax+1; ++j)
      V[i][j] = conserved_to_primitive(U[i][j]);

  applyMMSDirichletBCs(mmsCase, L, x_cell, y_cell);
}

void computeResidualMMS(
  int mmsCase,
  double gamma,
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell,
  const std::vector<std::vector<Primitive>>& V,
  std::vector<std::vector<Conserved>>& U,
  std::vector<std::vector<Conserved>>& R
) {
    int ni = imax + 2 * ghost;
    int nj = jmax + 2 * ghost;

    R.assign(ni, std::vector<Conserved>(nj, Conserved{0, 0, 0, 0}));

    for (int i = ghost; i <= imax + ghost; ++i) {
        for (int j = ghost; j < jmax + ghost; ++j) {
            Primitive VL = V[i-1][j];
            Primitive VR = V[i][j];
            double nx = nx_face_i[i][j];
            double ny = ny_face_i[i][j];
            double A = A_face_i[i][j];

            Conserved F = computeNumericalFlux(VL, VR, nx, ny, FluxScheme::VanLeer);
            R[i-1][j] -= F * A;
            R[i][j]   += F * A;
        }
    }

    for (int j = ghost; j <= jmax + ghost; ++j) {
        for (int i = ghost; i < imax + ghost; ++i) {
            Primitive VL = V[i][j-1];
            Primitive VR = V[i][j];
            double nx = nx_face_j[i][j];
            double ny = ny_face_j[i][j];
            double A = A_face_j[i][j];

            Conserved G = computeNumericalFlux(VL, VR, nx, ny, FluxScheme::VanLeer);
            R[i][j-1] -= G * A;
            R[i][j]   += G * A;
        }
    }
}




