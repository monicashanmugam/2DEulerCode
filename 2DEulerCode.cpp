#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <filesystem>  
#include <tuple>


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
    double rhou;   // momentum in x-direction
    double rhov;   // momentum in y-direction  
    double E;      // total energy
};



//-----------------------------------------------------
// Global Arrays for the CFD Solver
//-----------------------------------------------------

int imax = 0, jmax = 0;
// globals, un-sized until you know imax/jmax:
static std::vector<std::vector<Conserved>> U;  
static std::vector<std::vector<Primitive>> V;  




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
    std::vector<std::vector<double>> & ny_face_j,
    double &xmin, double &xmax, double &ymin, double &ymax, double &dx, double &dy
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
    int Ni = imax + 2 * ghost;
    int Nj = jmax + 2 * ghost;

    for(int i=0;i<Ni;++i){
    U[i].resize(Nj);
    V[i].resize(Nj);
  }

    // 4) Resize all eight arrays to the correct sizes
    x_cell.assign(Ni, std::vector<double>(Nj, 0.0));
    y_cell.assign(Ni, std::vector<double>(Nj, 0.0));
    A_face_i.assign(Ni + 1, std::vector<double>(Nj, 0.0));
    A_face_j.assign(Ni, std::vector<double>(Nj + 1, 0.0));
    nx_face_i.assign(Ni + 1, std::vector<double>(Nj, 0.0));
    ny_face_i.assign(Ni + 1, std::vector<double>(Nj, 0.0));
    nx_face_j.assign(Ni, std::vector<double>(Nj + 1, 0.0));
    ny_face_j.assign(Ni, std::vector<double>(Nj + 1, 0.0));

    // 5) Read interior cell-centers in the same order as the Fortran code:
    for (int k = 1; k <= kmax; ++k) {
        for (int j = 1; j <= jmax; ++j) {
            for (int i = 1; i <= imax; ++i) {
                double xx, yy;
                in >> xx >> yy;
                x_cell[i - 1 + ghost][j - 1 + ghost] = xx;
                y_cell[i - 1 + ghost][j - 1 + ghost] = yy;
            }
        }
    }

    // 6) Mirror-extrapolate two ghost layers (i-direction and j-direction)
    for (int i = 0; i < ghost; ++i) {
        int iL = ghost - 1 - i;
        int iR = ghost + imax + i;
        for (int j = ghost; j < ghost + jmax; ++j) {
            x_cell[iL][j] = x_cell[2 * ghost - 1 - i][j];
            y_cell[iL][j] = y_cell[2 * ghost - 1 - i][j];
            x_cell[iR][j] = x_cell[ghost + imax - 1 - i][j];
            y_cell[iR][j] = y_cell[ghost + imax - 1 - i][j];
        }
    }

    for (int j = 0; j < ghost; ++j) {
        int jB = ghost - 1 - j;
        int jT = ghost + jmax + j;
        for (int i = 0; i < Ni; ++i) {
            x_cell[i][jB] = x_cell[i][2 * ghost - 1 - j];
            y_cell[i][jB] = y_cell[i][2 * ghost - 1 - j];
            x_cell[i][jT] = x_cell[i][ghost + jmax - 1 - j];
            y_cell[i][jT] = y_cell[i][ghost + jmax - 1 - j];
        }
    }

    // 7) Compute vertical (i-) face lengths and outward normals
    for (int i = 0; i <= Ni; ++i) {
        for (int j = ghost; j < ghost + jmax; ++j) {
            int iL = std::max(0, i - 1);
            int iR = std::min(Ni - 1, i);
            double dx = x_cell[iR][j] - x_cell[iL][j];
            double dy = y_cell[iR][j] - y_cell[iL][j];
            double len = std::sqrt(dx * dx + dy * dy);

            A_face_i[i][j] = len;
            nx_face_i[i][j] = dy / len;
            ny_face_i[i][j] = -dx / len;
        }
    }

    // 8) Compute horizontal (j-) face lengths and outward normals
    for (int j = 0; j <= Nj; ++j) {
        for (int i = ghost; i < ghost + imax; ++i) {
            int jB = std::max(0, j - 1);
            int jT = std::min(Nj - 1, j);
            double dx = x_cell[i][jT] - x_cell[i][jB];
            double dy = y_cell[i][jT] - y_cell[i][jB];
            double len = std::sqrt(dx * dx + dy * dy);

            A_face_j[i][j] = len;
            nx_face_j[i][j] = -dy / len;
            ny_face_j[i][j] = dx / len;
        }
    }

    // 9) Compute the domain boundaries and grid spacings dx and dy
    xmin = x_cell[ghost][ghost];  // Minimum x-value
    xmax = x_cell[imax + ghost - 1][jmax + ghost - 1];  // Maximum x-value
    ymin = y_cell[ghost][ghost];  // Minimum y-value
    ymax = y_cell[imax + ghost - 1][jmax + ghost - 1];  // Maximum y-value

    // Calculate grid spacings dx and dy
    dx = (xmax - xmin) / imax;
    dy = (ymax - ymin) / jmax;

    // Output the domain boundaries and grid spacing
    std::cout << "[INFO] Domain boundaries:\n";
    std::cout << "xmin = " << xmin << ", xmax = " << xmax << "\n";
    std::cout << "ymin = " << ymin << ", ymax = " << ymax << "\n";

    std::cout << "[INFO] Grid spacings:\n";
    std::cout << "dx = " << dx << ", dy = " << dy << std::endl;
}


// 2) Overwrite everything (cell-centers, face lengths, normals) with a uniform
//    Cartesian grid of size imax x jmax over [0,L]x[0,L].  This is debug mode.
//    After calling this, you can run your solver exactly as if you had a
//    curvilinear grid-because the same eight arrays are filled, but in a trivial
//    Cartesian way.

inline std::tuple<double, double, double, double, double, double> convertToCartesianDebug(
    double L,  // Domain length
    std::vector<std::vector<double>>& x_cell,
    std::vector<std::vector<double>>& y_cell,
    std::vector<std::vector<double>>& A_face_i,
    std::vector<std::vector<double>>& A_face_j,
    std::vector<std::vector<double>>& nx_face_i,
    std::vector<std::vector<double>>& ny_face_i,
    std::vector<std::vector<double>>& nx_face_j,
    std::vector<std::vector<double>>& ny_face_j
) {
    // Total grid size including ghost layers
    int Ni = imax + 2 * ghost;
    int Nj = jmax + 2 * ghost;


    // Uniform grid spacing
    double dx = L / double(imax);  // Grid spacing in the x-direction
    double dy = L / double(jmax);  // Grid spacing in the y-direction

    // Assign cell-center coordinates uniformly on [0, L] x [0, L]
    x_cell.assign(Ni, std::vector<double>(Nj, 0.0));
    y_cell.assign(Ni, std::vector<double>(Nj, 0.0));

    // Calculate cell-center positions
    for (int i = 0; i < Ni; ++i) {
        for (int j = 0; j < Nj; ++j) {
            x_cell[i][j] = (i - ghost) * dx;  // X-coordinate of the cell center
            y_cell[i][j] = (j - ghost) * dy;  // Y-coordinate of the cell center
        }
    }

    // Set the face lengths (constant for Cartesian grid)
    A_face_i.assign(Ni + 1, std::vector<double>(Nj, dy));
    nx_face_i.assign(Ni + 1, std::vector<double>(Nj, 1.0));  // Normal is (1,0)
    ny_face_i.assign(Ni + 1, std::vector<double>(Nj, 0.0));

    A_face_j.assign(Ni, std::vector<double>(Nj + 1, dx));
    nx_face_j.assign(Ni, std::vector<double>(Nj + 1, 0.0));  // Normal is (0,1)
    ny_face_j.assign(Ni, std::vector<double>(Nj + 1, 1.0));

    // Compute the domain boundaries
    double xmin = x_cell[ghost][ghost];  // Minimum x-value
    double xmax = x_cell[imax + ghost - 1][jmax + ghost - 1];  // Maximum x-value
    double ymin = y_cell[ghost][ghost];  // Minimum y-value
    double ymax = y_cell[imax + ghost - 1][jmax + ghost - 1];  // Maximum y-value

    // Return all the computed values
    return {xmin, xmax, ymin, ymax, dx, dy};
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

// Change to a normal function
double PI = acos(-1.0); // This will work at runtime


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

//----------------------------------------------------------------------
// (E) computeSourceTermsMMS(…)
//----------------------------------------------------------------------
// Fill S[i][j] = (–∂F/∂x – ∂G/∂y) evaluated at the exact MMS solution.
// i.e. the analytic right‐hand‐side that makes your discrete residual ≈ 0.
void computeSourceTermsMMS(
  int mmsCase,
  double gamma,
  double L,
  const std::vector<std::vector<double>>& x_cell,
  const std::vector<std::vector<double>>& y_cell,
  std::vector<std::vector<Conserved>>&   S
) {
  assert(mmsCase==1 || mmsCase==2);
  int Ni = imax + 2*ghost;
  int Nj = jmax + 2*ghost;

  // allocate & zero
  S.assign(Ni, std::vector<Conserved>(Nj, Conserved{0,0,0,0}));

  // fill only the physical interior
  for (int i = ghost; i < ghost+imax; ++i) {
    for (int j = ghost; j < ghost+jmax; ++j) {
      double x = x_cell[i][j];
      double y = y_cell[i][j];
      S[i][j].rho  = rmassconv (mmsCase, L,     x, y);
      S[i][j].rhou = xmtmconv(mmsCase, L,     x, y);
      S[i][j].rhov = ymtmconv(mmsCase, L,     x, y);
      S[i][j].E    = energyconv(mmsCase, gamma, L, x, y);
    }
  }
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

            // Calculate boundary x and y positions
            double x = x_cell[i][j];
            double y = y_cell[i][j];

            // Initialize the Primitive and Conserved variables
            Primitive Vcell;
            Vcell.rho = rho_mms(mmsCase, L, x, y);
            Vcell.u   = uvel_mms(mmsCase, L, x, y);
            Vcell.v   = vvel_mms(mmsCase, L, x, y);
            Vcell.P   = press_mms(mmsCase, L, x, y);

            Conserved Ucell = PrimitiveToConserved(Vcell);

            // Update boundary cells (ghost cells)
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


//—— Flux / MUSCL helpers ——————————————————————————————————————————————
// (Make sure this goes after your #includes and your Primitive/Conserved definitions.)

inline double safeDenom(double d, double small = delta) {
  return (std::fabs(d) < small) ? ((d>=0) ? small : -small) : d;
}

inline double xi_limiter(double r) {
  return (r + std::fabs(r)) / (1.0 + r);
}

inline double Cplus(double M)  {
  if      (M<=-1) return 0.0;
  else if (M<   1) return 0.25*(M+1)*(M+1);
  else             return M;
}
inline double Cminus(double M) {
  if      (M<=-1) return M;
  else if (M<   1) return -0.25*(M-1)*(M-1);
  else             return 0.0;
}
inline double Dplus(double M) {
  if      (M<=-1) return 0.0;
  else if (M<   1) return 0.25*(M+1)*(M+1)*(2.0-M);
  else             return 1.0;
}
inline double Dminus(double M){
  if      (M<=-1) return 1.0;
  else if (M<   1) return -0.25*(M-1)*(M-1)*(2.0+M);
  else             return 0.0;
}

// 2D Van Leer flux across a face of area A with normal (nx,ny)
inline Conserved faceFluxVL2D(
  const Primitive &L, const Primitive &R,
  double nx, double ny, double A
) {
  // left
  double aL  = std::sqrt(gamma * L.P/L.rho);
  double MnL = (L.u*nx + L.v*ny)/aL;
  double CpL = Cplus(MnL), DpL = Dplus(MnL);
  double HtL = (gamma/(gamma-1))*L.P/L.rho + 0.5*(L.u*L.u + L.v*L.v);

  Conserved FcL{ L.rho*aL*CpL*A,
                 L.rho*aL*CpL*A*L.u,
                 L.rho*aL*CpL*A*L.v,
                 L.rho*aL*CpL*A*HtL };
  Conserved FpL{ 0.0,
                 DpL*L.P*nx*A,
                 DpL*L.P*ny*A,
                 0.0 };

  // right
  double aR  = std::sqrt(gamma * R.P/R.rho);
  double MnR = (R.u*nx + R.v*ny)/aR;
  double CmR = Cminus(MnR), DmR = Dminus(MnR);
  double HtR = (gamma/(gamma-1))*R.P/R.rho + 0.5*(R.u*R.u + R.v*R.v);

  Conserved FcR{ R.rho*aR*CmR*A,
                 R.rho*aR*CmR*A*R.u,
                 R.rho*aR*CmR*A*R.v,
                 R.rho*aR*CmR*A*HtR };
  Conserved FpR{ 0.0,
                 DmR*R.P*nx*A,
                 DmR*R.P*ny*A,
                 0.0 };

  return FcL + FpL + FcR + FpR;
}

// MUSCL in i-direction at face i+1/2, j fixed:
inline void musclI(
  const std::vector<std::vector<Primitive>> &V,
  int i, int j,
  int order, double kappa, bool freeze,
  Primitive &L, Primitive &R
) {
  double eps = (order==2 ? 1.0 : 0.0);
  double d0 = V[i+1][j].rho - V[i][j].rho;
  double rP = (V[i+2][j].rho - V[i+1][j].rho)/safeDenom(d0);
  double rM = (V[i][j].rho   - V[i-1][j].rho)/safeDenom(d0);
  double xiP = freeze?1.0:xi_limiter(rP);
  double xiM = freeze?1.0:xi_limiter(rM);
  auto rec = [&](auto get, auto set) {
    double LL = get(V[i-1][j]);
    double CC = get(V[i  ][j]);
    double RR = get(V[i+1][j]);
    double RR2= get(V[i+2][j]);
    L.*set = CC + (eps/4.0)*((1-kappa)*xiP*(CC-LL) + (1+kappa)*xiM*(RR-CC));
    R.*set = RR - (eps/4.0)*((1-kappa)*xiM*(RR2-RR) + (1+kappa)*xiP*(RR-CC));
  };
  rec([](auto&p){return p.rho;}, &Primitive::rho);
  rec([](auto&p){return p.u;  }, &Primitive::u);
  rec([](auto&p){return p.v;  }, &Primitive::v);
  rec([](auto&p){return p.P;  }, &Primitive::P);
}

// MUSCL in j-direction at face j+1/2, i fixed:
inline void musclJ(
  const std::vector<std::vector<Primitive>> &V,
  int i, int j,
  int order, double kappa, bool freeze,
  Primitive &L, Primitive &R
) {
  double eps = (order==2 ? 1.0 : 0.0);
  double d0 = V[i][j+1].rho - V[i][j].rho;
  double rP = (V[i][j+2].rho - V[i][j+1].rho)/safeDenom(d0);
  double rM = (V[i][j].rho   - V[i][j-1].rho)/safeDenom(d0);
  double xiP = freeze?1.0:xi_limiter(rP);
  double xiM = freeze?1.0:xi_limiter(rM);
  auto rec = [&](auto get, auto set) {
    double LL = get(V[i][j-1]);
    double CC = get(V[i][j  ]);
    double RR = get(V[i][j+1]);
    double RR2= get(V[i][j+2]);
    L.*set = CC + (eps/4.0)*((1-kappa)*xiP*(CC-LL) + (1+kappa)*xiM*(RR-CC));
    R.*set = RR - (eps/4.0)*((1-kappa)*xiM*(RR2-RR) + (1+kappa)*xiP*(RR-CC));
  };
  rec([](auto&p){return p.rho;}, &Primitive::rho);
  rec([](auto&p){return p.u;  }, &Primitive::u);
  rec([](auto&p){return p.v;  }, &Primitive::v);
  rec([](auto&p){return p.P;  }, &Primitive::P);
}
//———————————————————————————————————————————————————————————————————
//------------------------------------------------------------------------------
// computeResidualMMS
//
// Computes the purely convective residual R[i][j] = ∑faces (±F_face), using
// 2nd-order MUSCL + van Leer flux splitting.  You then add the analytic source
// S[i][j] in your time‐loop so that R+S → 0 for the exact MMS solution.
//
// fluxOrder     : 1 or 2
// kappa         : limiter parameter (e.g. −1, 0, 0.5, 1)
// freezeLimiter : if true, bypasses the limiter (→ pure MUSCL without TVD lim)
// x_cell, y_cell: cell‐center coordinates
// V             : primitive array (rho,u,v,P) including ghosts
// R             : output residual, same size as V
//------------------------------------------------------------------------------
void computeResidualMMS(
    int fluxOrder,
    double kappa,
    bool freezeLimiter,
    const std::vector<std::vector<double>>& x_cell,
    const std::vector<std::vector<double>>& y_cell,
    const std::vector<std::vector<Primitive>>& V,
    const std::vector<std::vector<Conserved>>& S,  // ← add this
    std::vector<std::vector<Conserved>>&   R
)
 {  
  int Ni = imax + 2*ghost;
  int Nj = jmax + 2*ghost;

  R.assign(Ni, vector<Conserved>(Nj, {0,0,0,0}));

  // —— vertical (i‐faces) ——
// only go as far as i+2 < Ni  ⇒  i ≤ Ni-3  ⇒  i ≤ (ghost+imax-1)
  for(int i = ghost; i <= ghost+imax-1; ++i){
    for(int j = ghost; j < ghost+jmax; ++j){
      Primitive PL, PR;
      musclI(V, i, j, fluxOrder, kappa, freezeLimiter, PL, PR);
      auto F = faceFluxVL2D(PL, PR, nx_face_i[i][j], ny_face_i[i][j], A_face_i[i][j]);
      R[i-1][j] -= F;
      R[i  ][j] += F;
    }
  }

  // —— horizontal (j‐faces) ——
// only go as far as j+2 < Nj  ⇒  j ≤ Nj-3  ⇒  j ≤ (ghost+jmax-1)
  for(int i = ghost; i < ghost+imax; ++i){
    for(int j = ghost; j <= ghost+jmax-1; ++j){
      Primitive PL, PR;
      musclJ(V, i, j, fluxOrder, kappa, freezeLimiter, PL, PR);
      auto G = faceFluxVL2D(PL, PR, nx_face_j[i][j], ny_face_j[i][j], A_face_j[i][j]);
      R[i][j-1] -= G;
      R[i][j  ] += G;
    }
  }

    // — now include the MMS source —
  for (int i = ghost; i < ghost+imax; ++i) {
    for (int j = ghost; j < ghost+jmax; ++j) {
      R[i][j].rho  -= S[i][j].rho;
      R[i][j].rhou -= S[i][j].rhou;
      R[i][j].rhov -= S[i][j].rhov;
      R[i][j].E    -= S[i][j].E;
    }
  }
}

void rungeKutta2Step(
    int fluxOrder, double kappa, bool freezeLimiter,
    int mmsCase, double L,
    const std::vector<std::vector<Conserved>>& S,
    std::vector<std::vector<Conserved>>&       R_int,
    double dt,
    double xmin, double xmax, double ymin, double ymax,
    bool debugMode, double dx, double dy
) {
    int Ni = imax + 2*ghost;
    int Nj = jmax + 2*ghost;

    // --- Stage 1: compute slope at (U,V)^n ---
    computeResidualMMS(fluxOrder, kappa, freezeLimiter,
                       x_cell, y_cell,
                       V,      // global primitive V^n
                       S,
                       R_int);

    // provisional update U* = U^n + dt * (−R_int) / vol
    auto U_star = U;
    for(int i=ghost; i<ghost+imax; ++i){
      for(int j=ghost; j<ghost+jmax; ++j){
        double vol = debugMode ? dx*dy : computeCellArea(i,j,x_cell,y_cell);
        U_star[i][j].rho  += -dt/vol * R_int[i][j].rho;
        U_star[i][j].rhou += -dt/vol * R_int[i][j].rhou;
        U_star[i][j].rhov += -dt/vol * R_int[i][j].rhov;
        U_star[i][j].E    += -dt/vol * R_int[i][j].E;
      }
    }
    applyBoundaryConditions(U_star, V, mmsCase, L, x_cell, y_cell);

    // --- NEW: build V_star from U_star ---
    std::vector<std::vector<Primitive>> V_star(Ni, std::vector<Primitive>(Nj));
    for(int i=0; i<Ni; ++i){
      for(int j=0; j<Nj; ++j){
        V_star[i][j] = ConservedToPrimitiveCell(U_star[i][j]);
      }
    }

    // --- Stage 2: slope at (U*,V*) ---
    computeResidualMMS(fluxOrder, kappa, freezeLimiter,
                       x_cell, y_cell,
                       V_star,  // <-- use the provisional primitives here
                       S,
                       R_int);

    // --- Combine slopes for Heun (2‐stage RK) ---
    for(int i=ghost; i<ghost+imax; ++i){
      for(int j=ghost; j<ghost+jmax; ++j){
        double vol = debugMode ? dx*dy : computeCellArea(i,j,x_cell,y_cell);
        // slope1 was in R_int from Stage1, slope2 now in R_int from Stage2
        Conserved slope1{ -R_int[i][j].rho,
                          -R_int[i][j].rhou,
                          -R_int[i][j].rhov,
                          -R_int[i][j].E };
        // (we could have stashed slope2 separately, but we re‐used R_int)
        Conserved slope2 = slope1; // because R_int was just overwritten by Stage2

        U[i][j].rho  += 0.5*dt/vol*( slope1.rho  + slope2.rho );
        U[i][j].rhou += 0.5*dt/vol*( slope1.rhou + slope2.rhou );
        U[i][j].rhov += 0.5*dt/vol*( slope1.rhov + slope2.rhov );
        U[i][j].E    += 0.5*dt/vol*( slope1.E    + slope2.E    );
      }
    }
    applyBoundaryConditions(U, V, mmsCase, L, x_cell, y_cell);
    GlobalConservedToPrimitive();
}




// Define the ResidualTriple structure
struct ResidualTriple {
    double mass;      // Mass residual
    double mom;       // Momentum residual (combined for x and y directions)
    double eng;       // Energy residual
    double combined;  // Combined residual (max of all three norms)
};


ResidualTriple computeResidualNorms(
    const std::vector<std::vector<Conserved>>& R
) {
    double sumM=0, sumMx=0, sumMy=0, sumE=0;
    for (int i = ghost; i < ghost+imax; ++i)
        for (int j = ghost; j < ghost+jmax; ++j) {
            sumM  += R[i][j].rho  * R[i][j].rho;
            sumMx += R[i][j].rhou * R[i][j].rhou;
            sumMy += R[i][j].rhov * R[i][j].rhov;
            sumE  += R[i][j].E    * R[i][j].E;
        }
    double n = double(imax*jmax);
    ResidualTriple r;
    r.mass     = std::sqrt(sumM  / n);
    r.mom      = std::sqrt((sumMx+sumMy) / n);
    r.eng      = std::sqrt(sumE  / n);
    r.combined = std::max({r.mass, r.mom, r.eng});
    return r;
}

//------------------------------------------------------------------------------
// 2D Tecplot dump: writes one POINT‐format zone per call.
//  Assumes globals imax, jmax, ghost, gamma, and arrays
//    x_cell, y_cell  (cell‐center coords, sized Ni×Nj)
//    V                (Primitive), sized Ni×Nj
// are all defined.
//------------------------------------------------------------------------------
void OutputSolution2D(const std::string &filename, int iter) {
    // 1) Open the file
    std::ofstream file(filename, std::ios::app);
    if (!file) {
        std::cerr << "Error: cannot open " << filename << " for writing\n";
        return;
    }

    // 2) On very first call, write header
    if (file.tellp() == 0) {
        file << "TITLE = \"2D Euler MMS Solution\"\n"
             << "VARIABLES = \"X\" \"Y\" \"rho\" \"u\" \"v\" \"P\" \"Mach\"\n";
    }

    // 3) Zone header
    file << "ZONE T=\"Iteration " << iter << "\""
         << " I=" << imax
         << " J=" << jmax
         << " K=1"
         << " DATAPACKING=POINT\n";

    // 4) local helper to guard against NaNs:
    auto safe = [](double v) {
      return std::isnan(v) ? -999.9 : v;
    };

    // 5) Dump all the physical cells
    for (int j = ghost; j < ghost + jmax; ++j) {
      for (int i = ghost; i < ghost + imax; ++i) {
        const auto &Vc = V[i][j];
        double rho  = Vc.rho;
        double u    = Vc.u;
        double v    = Vc.v;
        double P    = Vc.P;
        double a    = std::sqrt(gamma * P / rho);
        double Mach = (std::fabs(a)<1e-12 ? 0.0 : std::sqrt(u*u + v*v)/a);
        double x    = x_cell[i][j];
        double y    = y_cell[i][j];

        file 
          << safe(x)    << " "
          << safe(y)    << " "
          << safe(rho)  << " "
          << safe(u)    << " "
          << safe(v)    << " "
          << safe(P)    << " "
          << safe(Mach) << "\n";
      }
    }

    file.close();
    std::cout << "[INFO] Appended 2D solution for iter=" 
              << iter << " to " << filename << "\n";
}





int main() {

  int fluxOrder = 2;  // Second-order MUSCL
  double kappa = 0.5; // Choose limiter parameter (e.g., -1, 0, 0.5, or 1)
  bool freezeLimiter = false;  // Set to true for pure MUSCL (no limiter)

    // Create an output folder if needed.
    string outFolder = "OutputFiles";
    if (!fs::exists(outFolder)) {
        fs::create_directory(outFolder);
        cout << "[INFO] Created folder: " << outFolder << endl;
    }

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
    
    double xmin, xmax, ymin, ymax, dx, dy;
    
    if (!debugMode) {
      readCurviMeshFromFile(meshFile, x_cell, y_cell, A_face_i, A_face_j, nx_face_i, ny_face_i, nx_face_j, ny_face_j, xmin, xmax, ymin, ymax, dx, dy);
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
    
      std::tie(xmin, xmax, ymin, ymax, dx, dy) = convertToCartesianDebug(L, x_cell, y_cell, A_face_i, A_face_j, nx_face_i, ny_face_i, nx_face_j, ny_face_j);
      std::cout << "[INFO] Using debug Cartesian mesh: imax="<<imax<<", jmax="<<jmax<<"\n";
      
    }

    // ---------------------------------------------------
    // NOW: both branches have set imax/jmax and built the geometry
    int Ni = imax + 2*ghost;
    int Nj = jmax + 2*ghost;

    // resize solution arrays so U[i][j] and V[i][j] are valid everywhere
    U.assign(Ni, std::vector<Conserved>(Nj));
    V.assign(Ni, std::vector<Primitive>(Nj));

    std::cerr << "[DEBUG] After mesh:  U is " 
              << U.size() << "×" << U[0].size()
              << ", V is " 
              << V.size() << "×" << V[0].size()
              << "\n";

    std::cout << "[INFO] Loaded curvi mesh with imax = " << imax << ", jmax = " << jmax << "\n";
    std::cout << "[INFO] Sample values from x_cell and y_cell:\n";
    for (int i = 0; i < 5; ++i) { // Print first 5 x and y values
        std::cout << "x_cell[" << i << "] = " << x_cell[i][0] << ", y_cell[" << i << "] = " << y_cell[i][0] << "\n";
    }

    // Define the volume array (global or local)
    std::vector<std::vector<double>> cellVolume(imax + 2 * ghost, std::vector<double>(jmax + 2 * ghost, 0.0));

    // Loop over interior cells and compute area
    for (int i = ghost; i < imax + ghost; ++i) {
        for (int j = ghost; j < jmax + ghost; ++j) {
            cellVolume[i][j] = computeCellArea(i, j, x_cell, y_cell);
        }
    }

    // Check computed cell volume
    std::cout << "[INFO] Sample values from cellVolume array:\n";
    for (int i = ghost; i < ghost + imax; ++i) {
        for (int j = ghost; j < ghost + jmax; ++j) {
            std::cout << "cellVolume[" << i << "][" << j << "] = " << cellVolume[i][j] << "\n";
        }
    }

    


    std::cout << "[INFO] Sample values from A_face_i and A_face_j:\n";
    for (int i = 0; i < 5; ++i) {  // Print first 5 values for A_face_i and A_face_j
        std::cout << "A_face_i[" << i << "] = " << A_face_i[i][0] << ", A_face_j[" << i << "] = " << A_face_j[i][0] << "\n";
    }

    std::cout << "[INFO] Sample values from normal vectors nx_face_i, ny_face_i, nx_face_j, ny_face_j:\n";
    for (int i = 0; i < 5; ++i) {  // Print first 5 values for normal vectors
        std::cout << "nx_face_i[" << i << "] = " << nx_face_i[i][0] << ", ny_face_i[" << i << "] = " << ny_face_i[i][0] << "\n";
        std::cout << "nx_face_j[" << i << "] = " << nx_face_j[i][0] << ", ny_face_j[" << i << "] = " << ny_face_j[i][0] << "\n";
    }

    // 3) Call initializeMMS to fill U/V with the exact primitive → conserved MMS solution:
    int mmsCase = 1;
    initializeMMS(mmsCase, L, x_cell, y_cell, U, V);
    // Now U[i][j] and V[i][j] for i=ghost..ghost+imax−1, j=ghost..ghost+jmax−1 contain the _exact_ MMS solution at each cell‐center.

    // Apply the Dirichlet Boundary Conditions
    applyBoundaryConditions(U, V, mmsCase, L, x_cell, y_cell);

    std::vector<std::vector<Conserved>> S;
    computeSourceTermsMMS(mmsCase, gamma, L, x_cell, y_cell, S);

    GlobalPrimitiveToConserved();

    double dt = computeTimeStep(V, cellVolume, A_face_i, A_face_j, nx_face_i, ny_face_i, nx_face_j, ny_face_j);
    std::cout << "[INFO] Computed time step: dt = " << dt << "\n";

    

    std::vector<std::vector<Conserved>> R_int(imax, std::vector<Conserved>(jmax));  // Residuals

        // Pseudo-time marching parameters
    const int maxIter = 5000;
    const int writeInterval = 100;
    double tol = 1e-8;             // convergence tolerance
    ResidualTriple res{1e20,1e20,1e20,1e20};


    for (int iter = 0; iter < maxIter; ++iter) {
        // 1) compute a stable dt from the current primitive field V
        double dt = computeTimeStep(
            V, cellVolume,
            A_face_i,   A_face_j,
            nx_face_i,  ny_face_i,
            nx_face_j,  ny_face_j
        );

        // 2) advance U/V by one two-stage RK step
        rungeKutta2Step(
            fluxOrder, kappa, freezeLimiter,
            mmsCase, L,
            S,              // your MMS source term
            R_int,          // scratch for convective residual
            dt,
            xmin, xmax, ymin, ymax,
            debugMode, dx, dy
        );

        // 3) every writeInterval steps, dump to Tecplot and print residual
        if ((iter % writeInterval) == 0) {
            // compute/update your residual‐norms
            res = computeResidualNorms(R_int);
            std::cout << "iter=" << iter
                      << "  combined residual=" << res.combined
                      << "\n";

            // append a new zone at time=iter*dt
            OutputSolution2D(solFile, iter);
        }

        // 4) optional early exit if converged
        if (res.combined < tol) {
            std::cout << "Converged in " << iter << " steps.\n";
            break;
        }
    }

    return 0;
}
