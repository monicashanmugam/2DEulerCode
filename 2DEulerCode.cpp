#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <filesystem>  
#pragma once


using namespace std;

namespace fs = std::filesystem;

//-----------------------------------------------------
// Global constants and parameters
//-----------------------------------------------------
constexpr double Rgas      = 287.0;  // [J/(kg·K)]
constexpr double gamma = 1.4;
constexpr double L = 1.0; //domain length
static constexpr int ghost = 2;
const double CFL = 0.1;    // CFL number for time step control
const double epsM = 0.01;  // Minimum Mach number allowed
const double tolerance = 1e-12;
const double delta = 1e-6; // For computing the r

//-----------------------------------------------------
// Primitive and Conserved definitions
//-----------------------------------------------------

// Primitive variables: density (ρ), velocity (u,v), pressure (P)
struct Primitive {
    double rho;  // density
    double u;    // x‐velocity component
    double v;    // y‐velocity component 
    double P;    // pressure
};

// Conserved variables: mass, momentum, energy
struct Conserved {
    double rho;    // mass
    double rhou;   // momentum in x‐direction
    double rhov;   // momentum in y‐direction  
    double E;      // total energy
};


//-----------------------------------------------------
// Global Arrays for the CFD Solver
//-----------------------------------------------------

// Define the globals declared in variables.hpp
int imax;
int jmax;
std::vector<std::vector<Conserved>> U;
std::vector<std::vector<Primitive>> V;

// Inline operator overloads for “Conserved”
inline Conserved operator+(const Conserved &a, const Conserved &b) {
    return { a.rho  + b.rho,
             a.rhou + b.rhou,
             a.rhov + b.rhov,
             a.E    + b.E };
}

inline Conserved operator-(const Conserved &a, const Conserved &b) {
    return { a.rho  - b.rho,
             a.rhou - b.rhou,
             a.rhov - b.rhov,
             a.E    - b.E };
}

inline Conserved operator*(double s, const Conserved &b) {
    return { s * b.rho,
             s * b.rhou,
             s * b.rhov,
             s * b.E };
}

inline Conserved operator/(const Conserved &a, double s) {
    return { a.rho  / s,
             a.rhou / s,
             a.rhov / s,
             a.E    / s };
}

inline Conserved& operator+=(Conserved& a, const Conserved& b) {
    a.rho  += b.rho;
    a.rhou += b.rhou;
    a.rhov += b.rhov;
    a.E    += b.E;
    return a;
}

inline Conserved& operator-=(Conserved& a, const Conserved& b) {
    a.rho  -= b.rho;
    a.rhou -= b.rhou;
    a.rhov -= b.rhov;
    a.E    -= b.E;
    return a;
}


//-----------------------------------------------------
// Conversion Functions (cell-wise)
//-----------------------------------------------------

// Convert a cell's primitive variables (Vcell) to conserved variables (Ucell)
// Now, U = [rho, rho*u, rho*et]
Conserved PrimitiveToConserved(const Primitive &Vcell) {
    double kinetic_energy = 0.5 * Vcell.rho * (Vcell.u * Vcell.u + Vcell.v * Vcell.v);
    double E_vol = Vcell.P / (gamma - 1.0) + kinetic_energy;
    
    Conserved Ucell;
    Ucell.rho  = Vcell.rho;
    Ucell.rhou = Vcell.rho * Vcell.u;
    Ucell.rhov = Vcell.rho * Vcell.v;  
    Ucell.E    = E_vol;
    return Ucell;
}

Primitive ConservedToPrimitiveCell(const Conserved &Ucell) {
    Primitive Vcell;
    Vcell.rho = Ucell.rho;

    if (fabs(Ucell.rho) < 1e-12) {
        Vcell.u = 0.0;
        Vcell.v = 0.0;
    } else {
        Vcell.u = Ucell.rhou / Ucell.rho;
        Vcell.v = Ucell.rhov / Ucell.rho;  // new
    }

    double kinetic_energy = 0.5 * (Vcell.u * Vcell.u + Vcell.v * Vcell.v);
    Vcell.P = (gamma - 1.0) * (Ucell.E - Ucell.rho * kinetic_energy);

    if (Vcell.P < 1e-8)
        Vcell.P = 1e-8;

    return Vcell;
}


//-----------------------------------------------------
// Global Conversion Routines
//-----------------------------------------------------

void GlobalConservedToPrimitive() {
    int ni = imax + 2 * ghost;  // number of cells in i-direction (including ghosts)
    int nj = jmax + 2 * ghost;  // number of cells in j-direction (including ghosts)
    V.resize(ni);
    for (int i = 0; i < ni; i++) {
        V[i].resize(nj);
    }

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            V[i][j] = ConservedToPrimitiveCell(U[i][j]);
        }
    }
    // ApplyLimitsToPrimitive();
}

void GlobalPrimitiveToConserved() {
    int ni = imax + 2 * ghost;
    int nj = jmax + 2 * ghost;
    U.resize(ni);
    for (int i = 0; i < ni; i++) {
        U[i].resize(nj);
    }

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            U[i][j] = PrimitiveToConserved(V[i][j]);
        }
    }
    // ApplyLimitsToConserved();
}

// Define the eight geometry arrays:
std::vector<std::vector<double>> x_cell;
std::vector<std::vector<double>> y_cell;
std::vector<std::vector<double>> A_face_i;
std::vector<std::vector<double>> A_face_j;
std::vector<std::vector<double>> nx_face_i;
std::vector<std::vector<double>> ny_face_i;
std::vector<std::vector<double>> nx_face_j;
std::vector<std::vector<double>> ny_face_j;

struct ErrorData2D {
    double dx;     // grid spacing in i-direction
    double dy;     // grid spacing in j-direction
    double eP;     // L2 error for pressure
    double eRho;   // L2 error for density
    double eU;     // L2 error for u-velocity component
    double eV;     // L2 error for v-velocity component
};

// 1) Read a curvilinear mesh from a .grd file (and fill the eight arrays
//    listed above).  On return, the global imax/jmax are set from the file
//    header.  We also automatically create two layers of ghost cells by
//    mirror-extrapolating the interior points.

inline void readCurviMeshFromFile(
    const std::string & filepath,
    std::vector<std::vector<double>> & x_cell,
    std::vector<std::vector<double>> & y_cell,
    std::vector<std::vector<double>> & A_face_i,
    std::vector<std::vector<double>> & A_face_j,
    std::vector<std::vector<double>> & nx_face_i,
    std::vector<std::vector<double>> & ny_face_i,
    std::vector<std::vector<double>> & nx_face_j,
    std::vector<std::vector<double>> & ny_face_j
) {
    // 1) Open the file
    std::ifstream in(filepath);
    if (!in) {
        std::cerr << "[ERROR] Cannot open mesh file: " << filepath << "\n";
        std::exit(1);
    }

    // 2) Read nzones, imax, jmax, kmax into globals
    int nzones, kmax;
    in >> nzones >> imax >> jmax >> kmax;
    assert(nzones == 1 && kmax == 2);

    // 3) Compute total array dimensions (including two ghost layers on each side)
    int Ni = imax + 2*ghost;
    int Nj = jmax + 2*ghost;

    // 4) Resize all eight arrays to the correct sizes
    x_cell   .assign(Ni, std::vector<double>(Nj, 0.0));
    y_cell   .assign(Ni, std::vector<double>(Nj, 0.0));
    A_face_i .assign(Ni+1, std::vector<double>(Nj, 0.0));
    A_face_j .assign(Ni,   std::vector<double>(Nj+1, 0.0));
    nx_face_i.assign(Ni+1, std::vector<double>(Nj, 0.0));
    ny_face_i.assign(Ni+1, std::vector<double>(Nj, 0.0));
    nx_face_j.assign(Ni,   std::vector<double>(Nj+1, 0.0));
    ny_face_j.assign(Ni,   std::vector<double>(Nj+1, 0.0));

    // 5) Read interior cell-centers in the same order as the Fortran code:
    //      for k = 1..kmax
    //        for j = 1..jmax
    //          for i = 1..imax
    //            read x(i,j) -> store into x_cell[i-1+ghost][j-1+ghost]
    //    then repeat for y(i,j).
    for (int k = 1; k <= kmax; ++k) {
        for (int j = 1; j <= jmax; ++j) {
            for (int i = 1; i <= imax; ++i) {
                double xx;
                in >> xx;
                x_cell[i - 1 + ghost][j - 1 + ghost] = xx;
            }
        }
    }
    for (int k = 1; k <= kmax; ++k) {
        for (int j = 1; j <= jmax; ++j) {
            for (int i = 1; i <= imax; ++i) {
                double yy;
                in >> yy;
                y_cell[i - 1 + ghost][j - 1 + ghost] = yy;
            }
        }
    }

    // 6) Mirror-extrapolate two ghost layers in the i-direction (left/right)
    for (int i = 0; i < ghost; ++i) {
        int iL = ghost - 1 - i;
        int iR = ghost + imax + i;
        for (int j = ghost; j < ghost + jmax; ++j) {
            x_cell[iL][j] = x_cell[2*ghost - 1 - i][j];
            y_cell[iL][j] = y_cell[2*ghost - 1 - i][j];
            x_cell[iR][j] = x_cell[ghost + imax - 1 - i][j];
            y_cell[iR][j] = y_cell[ghost + imax - 1 - i][j];
        }
    }

    // 7) Mirror-extrapolate two ghost layers in the j-direction (bottom/top)
    for (int j = 0; j < ghost; ++j) {
        int jB = ghost - 1 - j;
        int jT = ghost + jmax + j;
        for (int i = 0; i < Ni; ++i) {
            x_cell[i][jB] = x_cell[i][2*ghost - 1 - j];
            y_cell[i][jB] = y_cell[i][2*ghost - 1 - j];
            x_cell[i][jT] = x_cell[i][ghost + jmax - 1 - j];
            y_cell[i][jT] = y_cell[i][ghost + jmax - 1 - j];
        }
    }

    // 8) Compute vertical (i-) face lengths and outward normals
    for (int i = 0; i <= Ni; ++i) {
        for (int j = ghost; j < ghost + jmax; ++j) {
            int iL = std::max(0,   i - 1);
            int iR = std::min(Ni-1, i);
            double dx = x_cell[iR][j] - x_cell[iL][j];
            double dy = y_cell[iR][j] - y_cell[iL][j];
            double len = std::sqrt(dx*dx + dy*dy);

            A_face_i[i][j]   = len;
            // Outward normal: pointing from right-cell toward left-cell
            nx_face_i[i][j] =  dy / len;
            ny_face_i[i][j] = -dx / len;
        }
    }

    // 9) Compute horizontal (j-) face lengths and outward normals
    for (int j = 0; j <= Nj; ++j) {
        for (int i = ghost; i < ghost + imax; ++i) {
            int jB = std::max(0,   j - 1);
            int jT = std::min(Nj-1, j);
            double dx = x_cell[i][jT] - x_cell[i][jB];
            double dy = y_cell[i][jT] - y_cell[i][jB];
            double len = std::sqrt(dx*dx + dy*dy);

            A_face_j[i][j]   = len;
            // Outward normal: pointing from top-cell toward bottom-cell
            nx_face_j[i][j] = -dy / len;
            ny_face_j[i][j] =  dx / len;
        }
    }
}

// 2) Overwrite everything (cell-centers, face lengths, normals) with a uniform
//    Cartesian grid of size imax x jmax over [0,L]x[0,L].  This is debug mode.
//    After calling this, you can run your solver exactly as if you had a
//    curvilinear grid-because the same eight arrays are filled, but in a trivial
//    Cartesian way.

inline void convertToCartesianDebug(
    double L,
    std::vector<std::vector<double>> & x_cell,
    std::vector<std::vector<double>> & y_cell,
    std::vector<std::vector<double>> & A_face_i,
    std::vector<std::vector<double>> & A_face_j,
    std::vector<std::vector<double>> & nx_face_i,
    std::vector<std::vector<double>> & ny_face_i,
    std::vector<std::vector<double>> & nx_face_j,
    std::vector<std::vector<double>> & ny_face_j
) {
    // We assume imax, jmax, and ghost are global constants or globals
    int Ni = imax + 2*ghost;
    int Nj = jmax + 2*ghost;

    // Uniform spacing
    double dx = L / double(imax);
    double dy = L / double(jmax);

    // 1) Assign cell-center coordinates uniformly on [0,L]x[0,L]
    x_cell.assign(Ni, std::vector<double>(Nj, 0.0));
    y_cell.assign(Ni, std::vector<double>(Nj, 0.0));
    for (int i = 0; i < Ni; ++i) {
        for (int j = 0; j < Nj; ++j) {
            x_cell[i][j] = (i - ghost + 0.5) * dx;
            y_cell[i][j] = (j - ghost + 0.5) * dy;
        }
    }

    // 2) Every vertical face (i-face) has the same length = dy,
    //    and outward normal = ( +1, 0 ).  We choose +1 so that
    //    face-i index i points into the cell at (i-1).
    A_face_i.assign(Ni+1, std::vector<double>(Nj, dy));
    nx_face_i.assign(Ni+1, std::vector<double>(Nj, +1.0));
    ny_face_i.assign(Ni+1, std::vector<double>(Nj,  0.0));

    // 3) Every horizontal face (j-face) has the same length = dx,
    //    and outward normal = ( 0, +1 ).  We choose +1 so that
    //    face-j index j points into the cell at (j-1).
    A_face_j.assign(Ni,   std::vector<double>(Nj+1, dx));
    nx_face_j.assign(Ni,   std::vector<double>(Nj+1,  0.0));
    ny_face_j.assign(Ni,   std::vector<double>(Nj+1, +1.0));
}

//----------------------------------------------------------------------
// (A) Define supersonic and subsonic MMS constants
//----------------------------------------------------------------------

// Constants for MMS field construction
struct MmsParams {
    double rho0,   rho_x,   rho_y,   a_rho_x, a_rho_y;
    double u0,     u_x,     u_y,     a_u_x,   a_u_y;
    double v0,     v_x,     v_y,     a_v_x,   a_v_y;
    double p0,     p_x,     p_y,     a_p_x,   a_p_y;
};

// Supersonic constants (MMS case 1)
constexpr MmsParams mmsSup = {
    1.0,   0.15,  -0.10,  1.0,    0.50,       // rho
    800.0, 50.0, -30.0,  1.5,    0.60,       // u
    800.0, -75.0, 40.0,  0.5,    (2.0 / 3.0),// v
    100000.0, 20000.0, 50000.0, 2.0, 1.0     // p
};

// Subsonic constants (MMS case 2)
constexpr MmsParams mmsSub = {
    1.0,   0.15,  -0.10,  1.0,   0.50,       // rho
    70.0,  5.0,   -7.0,  1.5,   0.60,       // u
    90.0, -15.0,   8.5,  0.5,   (2.0 / 3.0),// v
    100000.0, 20000.0, 50000.0, 2.0, 1.0     // p
};

// Shared constant for π
constexpr double PI = std::acos(-1.0);


double rho_mms(int mmsCase, double L, double x, double y) {
  assert(mmsCase == 1 || mmsCase == 2);
  const auto& C = (mmsCase == 1 ? mmsSup : mmsSub);

  double xi = x / L;
  double eta = y / L;
  double lin = C.rho0 + C.rho_x * xi + C.rho_y * eta;
  double sin_term = C.a_rho_x * std::sin(PI * xi)
                  + C.a_rho_y * std::sin(PI * eta);
  return lin + sin_term;
}

double uvel_mms(int mmsCase, double L, double x, double y) {
  assert(mmsCase == 1 || mmsCase == 2);
  const auto& C = (mmsCase == 1 ? mmsSup : mmsSub);

  double xi = x / L;
  double eta = y / L;
  double lin = C.u0 + C.u_x * xi + C.u_y * eta;
  double sin_term = C.a_u_x * std::sin(PI * xi)
                  + C.a_u_y * std::sin(PI * eta);
  return lin + sin_term;
}

double vvel_mms(int mmsCase, double L, double x, double y) {
  assert(mmsCase == 1 || mmsCase == 2);
  const auto& C = (mmsCase == 1 ? mmsSup : mmsSub);

  double xi = x / L;
  double eta = y / L;
  double lin = C.v0 + C.v_x * xi + C.v_y * eta;
  double sin_term = C.a_v_x * std::sin(PI * xi)
                  + C.a_v_y * std::sin(PI * eta);
  return lin + sin_term;
}

double press_mms(int mmsCase, double L, double x, double y) {
  assert(mmsCase == 1 || mmsCase == 2);
  const auto& C = (mmsCase == 1 ? mmsSup : mmsSub);

  double xi = x / L;
  double eta = y / L;
  double lin = C.p0 + C.p_x * xi + C.p_y * eta;
  double sin_term = C.a_p_x * std::sin(PI * xi)
                  + C.a_p_y * std::sin(PI * eta);
  return lin + sin_term;
}


//Initialize the cells with the exact solution from MMS

void initializeMMS(int mmsCase, double L,
                   const std::vector<std::vector<double>>& x_cell,
                   const std::vector<std::vector<double>>& y_cell,
                   std::vector<std::vector<Conserved>>& U,
                   std::vector<std::vector<Primitive>>& V) {
    
    int Ni = imax + 2 * ghost;
    int Nj = jmax + 2 * ghost;

    // Resize primitive array
    V.assign(Ni, std::vector<Primitive>(Nj));

    // Step 1: Fill Primitive Variables
    for (int i = 0; i < Ni; ++i) {
        for (int j = 0; j < Nj; ++j) {
            double x = x_cell[i][j];
            double y = y_cell[i][j];

            Primitive Vcell;
            Vcell.rho = rho_mms(mmsCase, L, x, y);
            Vcell.u   = uvel_mms(mmsCase, L, x, y);
            Vcell.v   = vvel_mms(mmsCase, L, x, y);
            Vcell.P   = press_mms(mmsCase, L, x, y);
            V[i][j] = Vcell;
        }
    }

    // Step 2: Use your pre-defined function to convert to Conserved
    GlobalPrimitiveToConserved();

    std::cout << "[INFO] Initialized MMS primitive + conserved.\n";
}

void applyBoundaryConditions(
    std::vector<std::vector<Conserved>> &U,
    std::vector<std::vector<Primitive>> &V,
    int mmsCase,
    double L,
    const std::vector<std::vector<double>> &x_cell,
    const std::vector<std::vector<double>> &y_cell
) {
    int ni = imax + 2 * ghost;
    int nj = jmax + 2 * ghost;

    for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
            if (i >= ghost && i < ghost + imax && j >= ghost && j < ghost + jmax) continue;

            double x = x_cell[i][j];
            double y = y_cell[i][j];

            Primitive Vcell;
            Vcell.rho = rho_mms(mmsCase, L, x, y);
            Vcell.u   = uvel_mms(mmsCase, L, x, y);
            Vcell.v   = vvel_mms(mmsCase, L, x, y);
            Vcell.P   = press_mms(mmsCase, L, x, y);

            Conserved Ucell = PrimitiveToConserved(Vcell);

            V[i][j] = Vcell;
            U[i][j] = Ucell;
        }
    }
}

double computeTimeStep(
    const std::vector<std::vector<Primitive>>& V,
    const std::vector<std::vector<double>>& cellVolume,
    const std::vector<std::vector<double>>& A_face_i,
    const std::vector<std::vector<double>>& A_face_j,
    const std::vector<std::vector<double>>& nx_face_i,
    const std::vector<std::vector<double>>& ny_face_i,
    const std::vector<std::vector<double>>& nx_face_j,
    const std::vector<std::vector<double>>& ny_face_j
) {
    double dtMin = 1e10;

    for (int i = ghost; i < imax + ghost; ++i) {
        for (int j = ghost; j < jmax + ghost; ++j) {
            const Primitive& Vcell = V[i][j];
            double a = std::sqrt(gamma * Vcell.P / Vcell.rho);
            double u = Vcell.u;
            double v = Vcell.v;

            // Eigenvalues across each face
            double lambda_iL = std::abs(u * nx_face_i[i][j] + v * ny_face_i[i][j]) + a;
            double lambda_iR = std::abs(u * nx_face_i[i+1][j] + v * ny_face_i[i+1][j]) + a;

            double lambda_jB = std::abs(u * nx_face_j[i][j] + v * ny_face_j[i][j]) + a;
            double lambda_jT = std::abs(u * nx_face_j[i][j+1] + v * ny_face_j[i][j+1]) + a;

            double areaSum =
                lambda_iL * A_face_i[i][j] +
                lambda_iR * A_face_i[i+1][j] +
                lambda_jB * A_face_j[i][j] +
                lambda_jT * A_face_j[i][j+1];

            double dt_cell = CFL * cellVolume[i][j] / areaSum;

            if (dt_cell < dtMin) dtMin = dt_cell;
        }
    }

    return dtMin;
}

// Compute area of cell (i,j) using diagonals AC and BD
double computeCellArea(
    int i, int j,
    const std::vector<std::vector<double>>& x_cell,
    const std::vector<std::vector<double>>& y_cell
) {
    // Corners (centered at 4 points)
    double xa = x_cell[i  ][j  ];
    double ya = y_cell[i  ][j  ];
    double xb = x_cell[i+1][j  ];
    double yb = y_cell[i+1][j  ];
    double xc = x_cell[i+1][j+1];
    double yc = y_cell[i+1][j+1];
    double xd = x_cell[i  ][j+1];
    double yd = y_cell[i  ][j+1];

    // Diagonals
    double ac_x = xc - xa;
    double ac_y = yc - ya;
    double bd_x = xb - xd;
    double bd_y = yb - yd;

    // Cross product magnitude
    double area = 0.5 * std::abs(ac_x * bd_y - ac_y * bd_x);
    return area;
}


// ---------------------------------------------------------------------
// ComputeFaceFlux:
// Compute the flux at a face given left and right primitive states,
// using central averaging. (Flux returned is per unit area.)
// ---------------------------------------------------------------------
Conserved ComputeFaceFlux(int iface, const Primitive &VL, const Primitive &VR) {
    (void)iface;
    Primitive Vface;
    Vface.rho = 0.5 * (VL.rho + VR.rho);
    Vface.u   = 0.5 * (VL.u   + VR.u);
    Vface.P   = 0.5 * (VL.P   + VR.P);

    Conserved flux;
    flux.rho  = Vface.rho * Vface.u;                              // Mass flux: ρ*u
    flux.rhou = Vface.rho * Vface.u * Vface.u + Vface.P;            // Momentum flux: ρ*u² + P
    flux.E    = (gamma_ / (gamma_ - 1.0)) * Vface.P * Vface.u 
                + Vface.rho * (0.5 * Vface.u * Vface.u * Vface.u);    // Energy flux
    return flux;
}


//------------------------------------------------------------------------------
// Helper function: safeDenom
// Returns d if |d| >= delta, else returns delta with the same sign as d.
//------------------------------------------------------------------------------
inline double safeDenom(double d, double delta_val = delta) {
    return (fabs(d) < delta_val) ? ((d >= 0.0) ? delta_val : -delta_val) : d;
}


//------------------------------------------------------------------------------
// Van Leer limiter: ξ(r) = (r + |r|) / (1 + r)
//------------------------------------------------------------------------------
inline double xi(double r) {
    return (r + fabs(r)) / (1.0 + r);
}

//------------------------------------------------------------------------------
// ******************* FLUX SPLITTING FUNCTIONS *******************
// Following van Leer (1982) style splitting.
// The following functions compute split flux coefficients and the physical flux.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Flux Splitting Functions (Updated Version)
//------------------------------------------------------------------------------

// SIGN(A, B): returns A with the sign of B.
inline double SIGN(double A, double B) {
    return (B >= 0.0 ? std::fabs(A) : -std::fabs(A));
}

// Flow direction sensors.
inline double alpha_plus(double M)  { 
    // If M > 0, return 1; if M < 0, return 0.
    if (M > 0.0)
        return 1.0;
    else
        return 0.0;
}

inline double alpha_minus(double M) {
    // If M < 0, return 1; if M > 0, return 0.
    if (M < 0.0)
        return 1.0;
    else
        return 0.0;
}

// Beta sensor: subsonic vs. supersonic.
inline double beta_sensor(double M) {
    // beta = -1 if |M| < 1, and 0 if |M| >= 1.
    if (std::fabs(M) < 1.0)
        return -1.0;
    else
        return 0.0;
}

// Split Mach coefficients (for convective flux):
// For Cplus: 0 if M <= -1, equals M+ if -1 < M < 1, equals M if M >= 1.
// M+ = 0.25*(M+1)^2.
inline double Cplus(double M) {
    if (M <= -1.0)
        return 0.0;
    else if (M < 1.0)
        return 0.25 * (M + 1.0) * (M + 1.0);
    else  // M >= 1.0
        return M;
}

// For Cminus: equals M if M <= -1, equals M^- if -1 < M < 1, 0 if M >= 1.
// M^- = -0.25*(M-1)^2.
inline double Cminus(double M) {
    if (M <= -1.0)
        return M;
    else if (M < 1.0)
        return -0.25 * (M - 1.0) * (M - 1.0);
    else  // M >= 1.0
        return 0.0;
}

// Pressure flux coefficients:
// Dplus: 0 if M <= -1, equals (M+)*(-M+2) if -1 < M < 1, equals 1 if M >= 1.
// Here, M+ is defined as above in Cplus.
inline double Dplus(double M) {
    if (M <= -1.0)
        return 0.0;
    else if (M < 1.0)
        return (0.25 * (M + 1.0) * (M + 1.0)) * (-M + 2.0);
    else  // M >= 1.0
        return 1.0;
}

// Dminus: 1 if M <= -1, equals (M^-)*(-M-2) if -1 < M < 1, equals 0 if M >= 1.
// Here, M^- is defined as in Cminus.
inline double Dminus(double M) {
    if (M <= -1.0)
        return 1.0;
    else if (M < 1.0)
        return (-0.25 * (M - 1.0) * (M - 1.0)) * (-M - 2.0);
    else  // M >= 1.0
        return 0.0;
}

//------------------------------------------------------------------------------
// First-order Van Leer Flux Function.
// (This function is unchanged; it computes the flux given left and right states.)
//------------------------------------------------------------------------------
Conserved ComputeFaceFluxVanLeerFirstOrder(const Primitive &VL, const Primitive &VR)
{
    //--- Left state calculations ---
    double aL = std::sqrt(gamma_ * VL.P / VL.rho);  // local speed of sound (left)
    double ML = VL.u / aL;                          // Mach number (left)
    double C_p = Cplus(ML);                         // split (convective) coefficient from left
    double D_p = Dplus(ML);                         // pressure split coefficient from left

    // Instead of computing total_enthalpy, we compute a combined energy term:
    double Ht_left_sum = (gamma_/(gamma_ - 1.0)) * (VL.P/VL.rho) + 0.5 * VL.u * VL.u;

    Conserved F_conv_left;
    F_conv_left.rho  = VL.rho * aL * C_p;
    F_conv_left.rhou = VL.rho * aL * C_p * VL.u;
    F_conv_left.E    = VL.rho * aL * C_p * Ht_left_sum;  // updated energy flux from left

    Conserved F_press_left;
    F_press_left.rho  = 0.0;
    F_press_left.rhou = D_p * VL.P;
    F_press_left.E    = 0.0;
    
    //--- Right state calculations ---
    double aR = std::sqrt(gamma_ * VR.P / VR.rho);  // local speed of sound (right)
    double MR = VR.u / aR;                          // Mach number (right)
    double C_m = Cminus(MR);                        // split coefficient from right
    double D_m = Dminus(MR);                        // pressure split coefficient from right

    double Ht_right_sum =  (gamma_/(gamma_ - 1.0)) * (VR.P/VR.rho) + 0.5 * VR.u * VR.u;

    Conserved F_conv_right;
    F_conv_right.rho  = VR.rho * aR * C_m;
    F_conv_right.rhou = VR.rho * aR * C_m * VR.u;
    F_conv_right.E    = VR.rho * aR * C_m * Ht_right_sum;  // updated energy flux from right

    Conserved F_press_right;
    F_press_right.rho  = 0.0;
    F_press_right.rhou = D_m * VR.P;
    F_press_right.E    = 0.0;
    
    //--- Total flux is the sum of convective and pressure contributions ---
    Conserved flux;
    flux.rho  = F_conv_left.rho  + F_conv_right.rho;
    flux.rhou = F_conv_left.rhou + F_conv_right.rhou + F_press_left.rhou + F_press_right.rhou;
    flux.E    = F_conv_left.E    + F_conv_right.E;
    
    return flux;
}


//------------------------------------------------------------------------------
// Compute van Leer flux at face i+1/2 with MUSCL + van Leer limiter
//------------------------------------------------------------------------------
Conserved ComputeFaceFluxVanLeer(
    const std::vector<Primitive>& V,
    int i,                // left cell index for face i+1/2
    int order,            // 1 or 2
    double kappa,
    bool freezeLimiter
) {
    double eps = (order == 2 ? 1.0 : 0.0);

    // Compute limiter factors for this face (i+1/2)
    // Denominator based on central jump
    double dU0 = V[i+1].u - V[i].u;
    double denom0 = safeDenom(dU0);
    // r+ at i+1/2 and r- at i+1/2
    double r_plus   = (V[i+2].u - V[i+1].u) / denom0;
    double r_minus  = (V[i].u   - V[i-1].u ) / denom0;
    double xi_p     = freezeLimiter ? 1.0 : xi(r_plus);
    double xi_m     = freezeLimiter ? 1.0 : xi(r_minus);

    // For VR: need limiter at face i+3/2
    double dU1 = V[i+2].u - V[i+1].u;
    double denom1 = safeDenom(dU1);
    double r_minus_ip1 = (V[i+1].u - V[i].u) / denom1;
    double xi_m_ip1    = freezeLimiter ? 1.0 : xi(r_minus_ip1);

    // Reconstruct left state at face i+1/2 (MUSCL with limiter)
    Primitive VL;
    VL.rho = V[i].rho + (eps/4.0) * (
        (1.0 - kappa) * xi_p      * (V[i].rho   - V[i-1].rho) +
        (1.0 + kappa) * xi_m      * (V[i+1].rho - V[i].rho  )
    );
    VL.u   = V[i].u   + (eps/4.0) * (
        (1.0 - kappa) * xi_p      * (V[i].u     - V[i-1].u  ) +
        (1.0 + kappa) * xi_m      * (V[i+1].u   - V[i].u    )
    );
    VL.P   = V[i].P   + (eps/4.0) * (
        (1.0 - kappa) * xi_p      * (V[i].P     - V[i-1].P  ) +
        (1.0 + kappa) * xi_m      * (V[i+1].P   - V[i].P    )
    );

    // Reconstruct right state at face i+1/2
    Primitive VR;
    VR.rho = V[i+1].rho - (eps/4.0) * (
        (1.0 - kappa) * xi_m_ip1 * (V[i+2].rho - V[i+1].rho) +
        (1.0 + kappa) * xi_p     * (V[i+1].rho - V[i].rho  )
    );
    VR.u   = V[i+1].u   - (eps/4.0) * (
        (1.0 - kappa) * xi_m_ip1 * (V[i+2].u   - V[i+1].u  ) +
        (1.0 + kappa) * xi_p     * (V[i+1].u   - V[i].u    )
    );
    VR.P   = V[i+1].P   - (eps/4.0) * (
        (1.0 - kappa) * xi_m_ip1 * (V[i+2].P   - V[i+1].P  ) +
        (1.0 + kappa) * xi_p     * (V[i+1].P   - V[i].P    )
    );

    // Compute and return the van Leer split flux
    return ComputeFaceFluxVanLeerFirstOrder(VL, VR);
}



//---------------------------------------------------------------------
// Flux function: computes physical flux vector f = [rho*u, rho*u^2 + P, (E+P)*u]
// using the definition E = P/(gamma-1) + 0.5*rho*u^2.
//---------------------------------------------------------------------
Conserved Flux(const Primitive &Vstate) {
    Conserved f;
    f.rho = Vstate.rho * Vstate.u;
    f.rhou = Vstate.rho * Vstate.u * Vstate.u + Vstate.P;
    double E = (gamma_ / (gamma_ - 1.0)) * Vstate.P * Vstate.u + Vstate.rho * (0.5 * Vstate.u * Vstate.u * Vstate.u);
    f.E = E;
    return f;
}


//---------------------------------------------------------------------
// Roe Averaging: Computes Roe-averaged quantities at the interface.
//---------------------------------------------------------------------
void RoeAverages(const Primitive &VL, const Primitive &VR,
    double &rhoRoe, double &uRoe, double &hRoe, double &aRoe)
{
rhoRoe = sqrt(VL.rho * VR.rho);
double sqrtL = sqrt(VL.rho);
double sqrtR = sqrt(VR.rho);
uRoe = (sqrtL * VL.u + sqrtR * VR.u) / (sqrtL + sqrtR);
// Compute specific enthalpy h = E/rho + P/rho, with E = P/(γ-1)+0.5*rho*u^2.
double hL = (gamma_/(gamma_ - 1.0)) * (VL.P / VL.rho) + 0.5 * VL.u * VL.u;
double hR = (gamma_/(gamma_ - 1.0)) * (VR.P / VR.rho) + 0.5 * VR.u * VR.u;
hRoe = (sqrtL * hL + sqrtR * hR) / (sqrtL + sqrtR);
aRoe = sqrt((gamma_ - 1.0) * (hRoe - 0.5 * uRoe * uRoe));
}

//---------------------------------------------------------------------
// RoeEigenstructure: Computes the eigenstructure and applies an entropy fix
// for Roe's flux formulation. The eigenvalues are modified as follows:
//   |λ_j|_mod = |λ_j|               if |λ_j| >= 2*epsilon_roe*aRoe,
//             = (λ_j^2)/(4*epsilon_roe*aRoe)+epsilon_roe*aRoe otherwise.
//---------------------------------------------------------------------
void RoeEigenstructure(double rhoRoe, double uRoe, double hRoe, double aRoe,
                       vector<vector<double>> &r, vector<double> &lambda)
{
    // Resize the eigenvector matrix and eigenvalue vector.
    r.resize(3, vector<double>(3, 0.0));
    lambda.resize(3, 0.0);
    
    // Compute right eigenvectors for the 1D Euler equations:
    // r1 = [1, uRoe, 0.5*uRoe^2]
    r[0][0] = 1.0;
    r[0][1] = uRoe;
    r[0][2] = 0.5 * uRoe * uRoe;
    // r2 = (rhoRoe/(2*aRoe))*[1, uRoe+aRoe, hRoe+uRoe*aRoe]
    double factor = rhoRoe / (2.0 * aRoe);
    r[1][0] = factor;
    r[1][1] = factor * (uRoe + aRoe);
    r[1][2] = factor * (hRoe + uRoe * aRoe);
    // r3 = - (rhoRoe/(2*aRoe))*[1, uRoe-aRoe, hRoe-uRoe*aRoe]
    r[2][0] = -factor;
    r[2][1] = -factor * (uRoe - aRoe);
    r[2][2] = -factor * (hRoe - uRoe * aRoe);
    
    // Define epsilon_roe (typically 0.1 or smaller).
    double epsilon_roe = 0.1;
    
    // Compute the raw eigenvalues.
    double lam0 = uRoe;         // corresponding to r1
    double lam1 = uRoe + aRoe;    // corresponding to r2
    double lam2 = uRoe - aRoe;    // corresponding to r3
    
    // A lambda fixer lambda_mod = (λ^2)/(4*epsilon_roe*aRoe) + epsilon_roe*aRoe if |λ| < 2*epsilon_roe*aRoe.
    auto fixLambda = [epsilon_roe, aRoe](double lam) -> double {
        double absLam = fabs(lam);
        if (absLam < 2 * epsilon_roe * aRoe) {
            return (lam * lam) / (4 * epsilon_roe * aRoe) + epsilon_roe * aRoe;
        } else {
            return absLam;
        }
    };
    
    // Store modified absolute eigenvalues.
    lambda[0] = fixLambda(lam0);
    lambda[1] = fixLambda(lam1);
    lambda[2] = fixLambda(lam2);
}


//---------------------------------------------------------------------
// Compute Delta Weights (Δw_j) for Roe flux formulation.
// Here we follow:
// deltaw1 = Δρ - ΔP/(aRoe^2)
// deltaw2 = Δu + ΔP/(rhoRoe*aRoe)
// deltaw3 = Δu - ΔP/(rhoRoe*aRoe)
// where Δu = u_R - u_L, Δρ = ρ_R - ρ_L, ΔP = P_R - P_L.
//---------------------------------------------------------------------
void ComputeDeltaWeights(const Primitive &VL, const Primitive &VR, double aRoe, double rhoRoe,
            double &dw1, double &dw2, double &dw3)
{
double delta_rho = VR.rho - VL.rho;
double deltaP = VR.P - VL.P;
double delta_u = VR.u - VL.u;
dw1 = delta_rho - (deltaP / (aRoe * aRoe));
dw2 = delta_u + (deltaP / (rhoRoe * aRoe));
dw3 = delta_u - (deltaP / (rhoRoe * aRoe));
}

//---------------------------------------------------------------------
// Compute Roe flux at face i+1/2 using MUSCL + Van Leer limiter
//---------------------------------------------------------------------
Conserved ComputeFaceFluxRoe(
    const std::vector<Primitive>& V,
    int i,                // left cell index for face i+1/2
    int order,            // 1 or 2
    double kappa,
    bool freezeLimiter
) {
    double eps = (order == 2 ? 1.0 : 0.0);

    // Compute limiter factors for this face (i+1/2)
    // Denominator based on central jump
    double dU0 = V[i+1].u - V[i].u;
    double denom0 = safeDenom(dU0);
    // r+ at i+1/2 and r- at i+1/2
    double r_plus   = (V[i+2].u - V[i+1].u) / denom0;
    double r_minus  = (V[i].u   - V[i-1].u ) / denom0;
    double xi_p     = freezeLimiter ? 1.0 : xi(r_plus);
    double xi_m     = freezeLimiter ? 1.0 : xi(r_minus);

    // For VR: need limiter at face i+3/2
    double dU1 = V[i+2].u - V[i+1].u;
    double denom1 = safeDenom(dU1);
    double r_minus_ip1 = (V[i+1].u - V[i].u) / denom1;
    double xi_m_ip1    = freezeLimiter ? 1.0 : xi(r_minus_ip1);

    // Reconstruct left state at face i+1/2 (MUSCL with limiter)
    Primitive VL;
    VL.rho = V[i].rho + (eps/4.0) * (
        (1.0 - kappa) * xi_p      * (V[i].rho   - V[i-1].rho) +
        (1.0 + kappa) * xi_m      * (V[i+1].rho - V[i].rho  )
    );
    VL.u   = V[i].u   + (eps/4.0) * (
        (1.0 - kappa) * xi_p      * (V[i].u     - V[i-1].u  ) +
        (1.0 + kappa) * xi_m      * (V[i+1].u   - V[i].u    )
    );
    VL.P   = V[i].P   + (eps/4.0) * (
        (1.0 - kappa) * xi_p      * (V[i].P     - V[i-1].P  ) +
        (1.0 + kappa) * xi_m      * (V[i+1].P   - V[i].P    )
    );

    // Reconstruct right state at face i+1/2
    Primitive VR;
    VR.rho = V[i+1].rho - (eps/4.0) * (
        (1.0 - kappa) * xi_m_ip1 * (V[i+2].rho - V[i+1].rho) +
        (1.0 + kappa) * xi_p     * (V[i+1].rho - V[i].rho  )
    );
    VR.u   = V[i+1].u   - (eps/4.0) * (
        (1.0 - kappa) * xi_m_ip1 * (V[i+2].u   - V[i+1].u  ) +
        (1.0 + kappa) * xi_p     * (V[i+1].u   - V[i].u    )
    );
    VR.P   = V[i+1].P   - (eps/4.0) * (
        (1.0 - kappa) * xi_m_ip1 * (V[i+2].P   - V[i+1].P  ) +
        (1.0 + kappa) * xi_p     * (V[i+1].P   - V[i].P    )
    );

    // Compute Roe averages
    double rhoRoe, uRoe, hRoe, aRoe;
    RoeAverages(VL, VR, rhoRoe, uRoe, hRoe, aRoe);

    // Compute the Roe eigenstructure and eigenvalues
    std::vector<std::vector<double>> r;
    std::vector<double> lambda;
    RoeEigenstructure(rhoRoe, uRoe, hRoe, aRoe, r, lambda);

    // Compute Delta Weights (deltaw1, deltaw2, deltaw3)
    double dw1, dw2, dw3;
    ComputeDeltaWeights(VL, VR, aRoe, rhoRoe, dw1, dw2, dw3);

    // Dissipation terms
    double diss_rho = std::fabs(lambda[0]) * dw1 * r[0][0] + 
                      std::fabs(lambda[1]) * dw2 * r[1][0] + 
                      std::fabs(lambda[2]) * dw3 * r[2][0];
    double diss_rhou = std::fabs(lambda[0]) * dw1 * r[0][1] + 
                        std::fabs(lambda[1]) * dw2 * r[1][1] + 
                        std::fabs(lambda[2]) * dw3 * r[2][1];
    double diss_E = std::fabs(lambda[0]) * dw1 * r[0][2] + 
                    std::fabs(lambda[1]) * dw2 * r[1][2] + 
                    std::fabs(lambda[2]) * dw3 * r[2][2];

    // Flux computation with dissipation
    Conserved flux;
    flux.rho  = 0.5 * (Flux(VL).rho + Flux(VR).rho) - 0.5 * diss_rho;
    flux.rhou = 0.5 * (Flux(VL).rhou + Flux(VR).rhou) - 0.5 * diss_rhou;
    flux.E    = 0.5 * (Flux(VL).E + Flux(VR).E) - 0.5 * diss_E;

    return flux;
}

// Define a structure to hold the three residual norms and a combined value.
struct ResidualTriple {
    double mass;
    double mom;
    double eng;
    double combined;  // For example, the maximum of the three.
};


//------------------------------------------------------------------------------
// static double initMass = 0.0, initMom = 0.0, initEng = 0.0;
// static ResidualTriple lastResidual = {0.0, 0.0, 0.0, 1e300};
// Global residual vector for interior cells (size = imax)
static vector<Conserved> R_int;

//------------------------------------------------------------------------------
// Initial residual norms and storage for lastResidual



//---------------------------------------------------------------------
// Updated UpdateSolution: Using van Leer upwind flux only (no JST dissipation)
// Returns a ResidualTriple structure with normalized residuals.
//---------------------------------------------------------------------
ResidualTriple UpdateSolution(double dt) {

    static int    iterCount      = 0;
    static bool   switchedTo2nd  = false;
    // // static double lastCheckRes   = 1e300;
    // // const int     holdIters      = 5000;   // initial first-order iterations
    // // const int     checkInterval  = 1000;   // check convergence frequency
    // // const double  tolDrop        = 1e-3;   // require 0.1% drop
    const int     maxIterSwitch  = 20000;  // Force switch to second-order after this many iterations

    // ++iterCount;



    // // Check for convergence rate and decide whether to switch to second-order flux
    // if (!switchedTo2nd && iterCount >= holdIters && (iterCount % checkInterval) == 0) {
    //     double prevRes = lastCheckRes;
    //     double currRes = lastResidual.combined;  // assumed stored

    //     // Check if the residual is not dropping fast enough
    //     if (currRes > prevRes * (1.0 - tolDrop)) {
    //         switchedTo2nd = true;
    //         std::cout << "[INFO] Switching to 2nd-order flux at iteration "
    //                   << iterCount << ".\n";
    //     }
    //     lastCheckRes = currRes;
    // }

    // Force the switch to second-order after maxIterSwitch iterations
    if (iterCount >= maxIterSwitch && !switchedTo2nd) {
        switchedTo2nd = true;
        std::cout << "[INFO] Forcing switch to 2nd-order flux at iteration "
                  << iterCount << ".\n";
    }

    // Check if we should freeze the limiter (and if not already frozen)
    // Inside UpdateSolution, after incrementing the iteration counter:
    if (iterationCount >= 50000 && !freezeLimiter) {
        freezeLimiter = true;
        cout << "[INFO] Limiters frozen: now updating solution without applying limiter functions." << endl;
    }

    // Determine fluxOrder: start 1st-order, then switch based on convergence
    int fluxOrder = switchedTo2nd ? 2 : 1;

    // int fluxOrder = 1;
    
    // Set the flux order and kappa.
    double kappa = -1;  // Choose appropriate value: e.g., -1, 0, 0.5, or 1
    int numPhysFaces = imax + 1;          // physical faces (indices 0 to imax)
    // Build the physical face flux vector F_int using van Leer upwind flux.
    vector<Conserved> F_int(numPhysFaces);
    // For face 0: interface between cell (ghost-1) and cell (ghost).
    // (Adjust indices if needed; here we assume V[ghost-1] is valid in the ghost region.)
    // F_int[0] = ComputeFaceFluxVanLeer(V, ghost - 1, fluxOrder, kappa);
    for (int j = 0; j < numPhysFaces; j++) {
        // When j=1, this is the face between V[ghost+0] and V[ghost+1],
        // so we pass ghost+0 as the left index.
        F_int[j] = ComputeFaceFluxRoe(V, ghost + j - 1, fluxOrder, kappa, freezeLimiter);
    }
    // For the last face (physical face index = numPhysFaces-1)
    // F_int[numPhysFaces - 1] = ComputeFaceFluxVanLeer(V, lastInterior, fluxOrder, kappa);
    
    
    // Now update the solution solely based on the difference in upwind fluxes.
    double dx = (xmax - xmin) / imax;
    vector<Conserved> U_new = U;
    R_int.resize(imax);
    
    // Loop over each interior (physical) cell.
    for (int j = 0; j < imax; j++) {
        int i = ghost + j;
        
        // Get left and right face fluxes from F_int.
        // Multiply each by the corresponding face area (A_face_phys).
        Conserved F_left = F_int[j];
        F_left.rho  *= A_face_phys[j];
        F_left.rhou *= A_face_phys[j];
        F_left.E    *= A_face_phys[j];
        
        Conserved F_right = F_int[j + 1];
        F_right.rho  *= A_face_phys[j + 1];
        F_right.rhou *= A_face_phys[j + 1];
        F_right.E    *= A_face_phys[j + 1];
        
        // Net flux leaving the cell is difference between right and left flux.
        Conserved netFlux = F_right - F_left;
        
        // Geometric source term (for quasi-1D nozzle)
        double areaDiff = A_face_phys[j + 1] - A_face_phys[j];
        double S_momentum = V[i].P * areaDiff / dx;
        Conserved S = {0.0, S_momentum, 0.0};
        
        netFlux.rhou -= S.rhou * dx;
        
        // Compute the cell volume (area*dx), where cell area is average of face areas.
        double cellArea   = 0.5 * (A_face_phys[j] + A_face_phys[j + 1]);
        double cellVolume = cellArea * dx;
        
        // Euler explicit update.
        U_new[i].rho  = U[i].rho  - (dt / cellVolume) * netFlux.rho;
        U_new[i].rhou = U[i].rhou - (dt / cellVolume) * netFlux.rhou;
        U_new[i].E    = U[i].E    - (dt / cellVolume) * netFlux.E;
        
        // Store the residual for the current cell.
        R_int[j] = netFlux;
    }
    
    // Update interior cells.
    for (int j = 0; j < imax; j++) {
        int i = ghost + j;
        U[i] = U_new[i];
    }
    
    // Update ghost cells (if needed) and convert new U to primitive variables.
    ApplyLimitsToConserved();
    GlobalConservedToPrimitive();
    
    // -------------------------------
    // Compute residual norms.
    // -------------------------------
    double sumMass = 0.0, sumMom = 0.0, sumEng = 0.0;
    for (int j = 0; j < imax; j++) {
        sumMass += R_int[j].rho * R_int[j].rho;
        sumMom  += R_int[j].rhou * R_int[j].rhou;
        sumEng  += R_int[j].E * R_int[j].E;
    }
    double currentMass = sqrt(sumMass / imax);
    double currentMom  = sqrt(sumMom / imax);
    double currentEng  = sqrt(sumEng / imax);
    
    static bool normBaseSet = false;
    static double initMass = 0.0, initMom = 0.0, initEng = 0.0;
    if (!normBaseSet) {
        initMass = currentMass;
        initMom  = currentMom;
        initEng  = currentEng;
        normBaseSet = true;
    }
    
    double resMass_norm = (fabs(initMass) < 1e-12 ? currentMass : currentMass / initMass);
    double resMom_norm  = (fabs(initMom)  < 1e-12 ? currentMom  : currentMom / initMom);
    double resEng_norm  = (fabs(initEng)  < 1e-12 ? currentEng  : currentEng / initEng);
    
    cout << "Residual Norms (Normalized): Mass = " << resMass_norm
         << ", Momentum = " << resMom_norm
         << ", Energy = " << resEng_norm << endl;
    
    double combinedResidual = std::max({resMass_norm, resMom_norm, resEng_norm});
    cout << "[INFO] Combined residual: " << combinedResidual << endl;
    
    ResidualTriple residualData;
    residualData.mass = resMass_norm;
    residualData.mom = resMom_norm;
    residualData.eng = resEng_norm;
    residualData.combined = combinedResidual;
    
    return residualData;
}



// Updated OutputSolution: Always writes a Tecplot zone header
// and then the interior solution data.
// We assume that global physical arrays x_cell_phys (size = imax)
// and A_face_phys (size = imax+1) have been set in SetGeometry,
// and that the full solution (in V) is stored for indices ghost ... ghost+imax-1.
void OutputSolution(const string &filename, int iter, bool cond) {
    // Open the file in append mode
    ofstream file(filename, ios::app);
    if (!file) {
        cerr << "Error: cannot open file " << filename << " for writing." << endl;
        return;
    }
    
    // If the file is empty, write the global header (TITLE and VARIABLES).
    // (This will only run on the very first call when the file is new.)
    if (file.tellp() == 0) {
        file << "TITLE = \"Quasi-1D Nozzle Solutions\"\n";
        file << "VARIABLES = \"x\" \"A_cell\" \"rho\" \"u\" \"P\" \"Mach\"\n";
    }
    
    // Write a zone header for this iteration.
    file << "ZONE T=\"Iteration " << iter << "\", I=" << imax << ", F=POINT\n";
    
    // Loop over interior (physical) cells. 
    // x_cell_phys holds cell-centers (length imax), and A_face_phys holds face-areas (length imax+1).
    for (int i = 0; i < imax; i++) {
        // Compute cell area as the average of the left and right face areas.
        double cellArea = 0.5 * (A_face_phys[i] + A_face_phys[i+1]);
        // Map physical cell index i to full domain index (cells start at index 'ghost').
        int fullIndex = ghost + i;
        double a = sqrt(gamma_ * V[fullIndex].P / V[fullIndex].rho);
        double Mach = (fabs(a) < 1e-12 ? 0.0 : V[fullIndex].u / a);
        // Output: cell center, cell area, density, velocity, pressure, Mach.
        file << x_cell_phys[i] << " " << cellArea << " " 
             << V[fullIndex].rho << " " << V[fullIndex].u << " " 
             << V[fullIndex].P << " " << Mach << "\n";
    }
    
    file.close();
    cout << "[INFO] Appended solution for iteration " << iter << " to " << filename << endl;
}


int main() {


    // Create an output folder if needed.
    string outFolder = "OutputFiles";
    if (!fs::exists(outFolder)) {
        fs::create_directory(outFolder);
        cout << "[INFO] Created folder: " << outFolder << endl;
    }

    // Open a file for error norms.
    ofstream tecFile("ErrorNormsTecplot.dat", ios::out);
    if (!tecFile) {
        cerr << "[ERROR] Could not open ErrorNormsTecplot.dat for writing!" << endl;
        return 1;
    }
    tecFile << "TITLE = \"Error Norms for 2D MMS Case (First Order)\"\n";
    tecFile << "VARIABLES = \"dx\" \"Pressure_Error\" \"Density_Error\" \"Velocity_Error\"\n";

    // Open (or create) a Tecplot file for solution output.
    // We want all zones (initial condition and every 100 iterations) in one file.
    std::string solFile = outFolder + "/MMS_Solution.dat";
    // Erase any existing file by opening in output mode once.
    std::ofstream solOut(solFile, std::ios::out);
    solOut.close();

    
    // 1) Construct the folder / file you want to read:
    const std::string meshFile = R"(C:\Users\monicashanmugam\OneDrive - Virginia Tech\Desktop\Virginia Tech\CFD\Project\Project_Files\Project_Files\Grids\curviliniear-grids\curv2d9.grd)";

    // 2) Decide whether you want to run in DEBUG (Cartesian) mode
    bool debugMode = true;  // set true if you want a simple Cartesian mesh   
    
    if (!debugMode) {
      // read full curvilinear geometry from file:
      readCurviMeshFromFile(meshFile, x_cell, y_cell, A_face_i, A_face_j, nx_face_i, ny_face_i, nx_face_j, ny_face_j);
      std::cout << "[INFO] Loaded curvi mesh with imax = " << imax << ", jmax = " << jmax << "\n";
    }
    else {
      // in debug mode, only read imax/jmax and then overwrite with Cartesian:
      {
        std::ifstream in(meshFile);
        int nz, kmax;
        in >> nz >> imax >> jmax >> kmax;
        assert(nz == 1 && kmax == 2);
      }
    
      convertToCartesianDebug( L, x_cell, y_cell, A_face_i, A_face_j, nx_face_i, ny_face_i, nx_face_j, ny_face_j );
      std::cout << "[INFO] Using debug Cartesian mesh: imax="<<imax<<", jmax="<<jmax<<"\n";
    }

    // Define the volume array (global or local)
    std::vector<std::vector<double>> cellVolume(imax + 2 * ghost, std::vector<double>(jmax + 2 * ghost, 0.0));

    // Loop over interior cells
    for (int i = ghost; i < imax + ghost; ++i) {
        for (int j = ghost; j < jmax + ghost; ++j) {
            cellVolume[i][j] = computeCellArea(i, j, x_cell, y_cell);
        }
    }


    // 3) Call initializeMMS to fill U/V with the exact primitive → conserved MMS solution:
    int mmsCase = 1;
    initializeMMS(mmsCase, L, x_cell, y_cell, U, V);

    // Now U[i][j] and V[i][j] for i=ghost..ghost+imax−1, j=ghost..ghost+jmax−1
      // contain the _exact_ MMS solution at each cell‐center.


    // Apply the Dirichlet Boundary Conditions
    applyBoundaryConditions(U, V, mmsCase, L, x_cell, y_cell);

    double dt = computeTimeStep(V, cellVolume, A_face_i, A_face_j, nx_face_i, ny_face_i, nx_face_j, ny_face_j);
    std::cout << "[INFO] Computed time step: dt = " << dt << std::endl;






      // --- Output the initial solution (iteration 0) ---
      OutputSolution(solFile, 0, true);

        // --- Run solver iterations ---
      int nmax = 120000;  // You can adjust the maximum iterations.
      vector<ResidualTriple> residualHistory;
      for (int niter = 0; niter < nmax; niter++) {
          double dt = ComputeTimeStep();
          ResidualTriple resData = UpdateSolution(dt);
          GlobalConservedToPrimitive();
          ApplyBoundaryConditions();
          cout << "[INFO] Iteration " << niter 
               << ": Residual (Mass, Mom, Eng) = (" 
               << resData.mass << ", " << resData.mom << ", " 
               << resData.eng << ")" << endl;
          residualHistory.push_back(resData);
          if (resData.mass < 1e-9 && resData.mom < 1e-9 && resData.eng < 1e-9) {
              cout << "[INFO] Convergence achieved at iteration " << niter 
                   << " with residuals: (" << resData.mass << ", " 
                   << resData.mom << ", " << resData.eng << ")" << endl;
              // Output the solution at convergence.
              OutputSolution(solFile, niter, true);
              break;
          }
          // Every 100 iterations, output the FVM solution as a new Tecplot zone.
          if (niter % 10000 == 0) {
            OutputSolution(solFile, niter, true);
          }
      }

      // --- Write residual history to a CSV file ---
      string residualFile = outFolder + "/Residuals_imax" + to_string(imax) + ".csv";
      ofstream resFile(residualFile, ios::out);
      if (!resFile) {
          cerr << "[ERROR] Cannot open file " << residualFile << " for writing." << endl;
      } else {
          resFile << "Iteration,Mass,Momentum,Energy\n";
          for (size_t i = 0; i < residualHistory.size(); i++) {
              resFile << i << "," << residualHistory[i].mass << "," 
                      << residualHistory[i].mom << "," << residualHistory[i].eng << "\n";
          }
          resFile.close();
          cout << "[INFO] Residual history written to " << residualFile << endl;
      }

      // --- Compute L2 error norms using cell-averaged exact solution ---
      double error_rho = 0.0, error_u = 0.0, error_P = 0.0;
      int N = imax; // number of physical cells.
      for (int i = 0; i < N; i++) {
          int fullIndex = ghost + i;
          double rho_h = V[fullIndex].rho;
          double u_h   = V[fullIndex].u;
          double P_h   = V[fullIndex].P;
          Primitive exactAvg = ExactSolutionAverageOverCell(i);
          double diff_rho = rho_h - exactAvg.rho;
          double diff_u   = u_h - exactAvg.u;
          double diff_P   = P_h - exactAvg.P;
          error_rho += diff_rho * diff_rho;
          error_u   += diff_u * diff_u;
          error_P   += diff_P * diff_P;
      }
      error_rho = sqrt(error_rho / N);
      error_u   = sqrt(error_u / N);
      error_P   = sqrt(error_P / N);

      // --- Compute dx and output the Tecplot zone for error norms ---
      double dx_val = (xmax - xmin) / imax;
      tecFile << "ZONE T=\"Mesh " << imax << "\", I=1, F=POINT\n";
      tecFile << dx_val << " " << error_P << " " << error_rho << " " << error_u << "\n";

      cout << "[INFO] Completed run for imax = " << imax 
             << " => dx = " << dx_val 
             << ", eP = " << error_P 
             << ", eRho = " << error_rho 
             << ", eU = " << error_u << endl;

           

      tecFile.close();
      cout << "[INFO] Wrote Tecplot file: ErrorNormsTecplot.dat" << endl;

      return 0;
}
