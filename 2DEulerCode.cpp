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

std::vector<std::vector<double>> x_node;
std::vector<std::vector<double>> y_node;


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
int kmax;
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


//------------------------------------------------------------------------------
// Global Limits for Primitive variables
//------------------------------------------------------------------------------
const double RHO_MIN = 1e-3;
const double RHO_MAX = 5.0;

const double U_MIN   = -1000.0;
const double U_MAX   = 2000.0;
const double V_MIN   = -1000.0;
const double V_MAX   = 2000.0;

const double P_MIN   = 1e-8;
const double P_MAX   = 1e6;

//------------------------------------------------------------------------------
// Global Limits for Conserved variables
//------------------------------------------------------------------------------
const double U1_MIN = RHO_MIN * 0.5;
const double U1_MAX = RHO_MAX * 2.0;

const double U2_MIN = RHO_MIN * U_MIN * 2.0;   // rho*u
const double U2_MAX = RHO_MAX * U_MAX * 2.0;

const double U3_MIN = RHO_MIN * V_MIN * 2.0;   // rho*v
const double U3_MAX = RHO_MAX * V_MAX * 2.0;

const double U4_MIN = P_MIN / (gamma - 1.0);  // energy
const double U4_MAX = P_MAX / (gamma - 1.0)
                    + 0.5 * RHO_MAX * (U_MAX*U_MAX + V_MAX*V_MAX);


// call this after you fill V (primitive) but before you convert to U:
void ApplyLimitsToPrimitive(const std::string &when,
                            std::vector<std::vector<Primitive>> &V) {
    int ni = V.size();
    if (ni==0) return;
    int nj = V[0].size();

    // store up to 5 clamped locations
    std::vector<std::pair<int,int>> clampedCells;
    clampedCells.reserve(5);

    for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
            bool clamped = false;
            auto &cell = V[i][j];

            double new_rho = std::clamp(cell.rho, RHO_MIN, RHO_MAX);
            if (new_rho != cell.rho) { clamped = true; cell.rho = new_rho; }

            double new_u = std::clamp(cell.u, U_MIN, U_MAX);
            if (new_u != cell.u) { clamped = true; cell.u = new_u; }

            double new_v = std::clamp(cell.v, V_MIN, V_MAX);
            if (new_v != cell.v) { clamped = true; cell.v = new_v; }

            double new_P = std::clamp(cell.P, P_MIN, P_MAX);
            if (new_P != cell.P) { clamped = true; cell.P = new_P; }

            if (clamped && clampedCells.size() < 5) {
                clampedCells.emplace_back(i,j);
            }
        }
    }

    if (!clampedCells.empty()) {
        std::cerr << "[WARNING] Clamped primitive vars"
                  << (when.empty() ? "" : " ("+when+")")
                  << " in " 
                  << /* total count */ /* you can track a counter if you want fully accurate count */
                  clampedCells.size()
                  << " cell"
                  << (clampedCells.size()>1?"s":"")
                  << ", e.g. ";
        for (auto [i,j] : clampedCells) {
            std::cerr << "("<<i<<","<<j<<") ";
        }
        if (clampedCells.size()==5) std::cerr << "…";
        std::cerr << "\n";
    }
}

// same idea for conserved U:
void ApplyLimitsToConserved(const std::string &when,
                            std::vector<std::vector<Conserved>> &U) {
    int ni = U.size();
    if (ni==0) return;
    int nj = U[0].size();

    std::vector<std::pair<int,int>> clampedCells;
    clampedCells.reserve(5);

    for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
            bool clamped = false;
            auto &cell = U[i][j];

            // define these U?_MIN/U?_MAX appropriately…
            double new_rho  = std::clamp(cell.rho,  U1_MIN, U1_MAX);
            if (new_rho != cell.rho) { clamped = true; cell.rho = new_rho; }

            double new_rhou = std::clamp(cell.rhou, U2_MIN, U2_MAX);
            if (new_rhou != cell.rhou) { clamped = true; cell.rhou = new_rhou; }

            double new_rhov = std::clamp(cell.rhov, U3_MIN, U3_MAX);
            if (new_rhov != cell.rhov) { clamped = true; cell.rhov = new_rhov; }

            double new_E    = std::clamp(cell.E,    U4_MIN, U4_MAX);
            if (new_E != cell.E) { clamped = true; cell.E = new_E; }

            if (clamped && clampedCells.size() < 5) {
                clampedCells.emplace_back(i,j);
            }
        }
    }

    if (!clampedCells.empty()) {
        std::cerr << "[WARNING] Clamped conserved vars"
                  << (when.empty() ? "" : " ("+when+")")
                  << " in " 
                  << clampedCells.size()
                  << " cell"
                  << (clampedCells.size()>1?"s":"")
                  << ", e.g. ";
        for (auto [i,j] : clampedCells) {
            std::cerr << "("<<i<<","<<j<<") ";
        }
        if (clampedCells.size()==5) std::cerr << "…";
        std::cerr << "\n";
    }
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
    int ni = imax + 2*ghost, nj = jmax + 2*ghost;
    V.resize(ni);
    for(int i=0;i<ni;++i) V[i].resize(nj);

    for(int i=0;i<ni;++i)
      for(int j=0;j<nj;++j)
        V[i][j] = ConservedToPrimitiveCell(U[i][j]);

    ApplyLimitsToPrimitive("after U→V conversion", V);
}

void GlobalPrimitiveToConserved() {
    int ni = imax + 2*ghost, nj = jmax + 2*ghost;
    U.resize(ni);
    for(int i=0;i<ni;++i) U[i].resize(nj);

    for(int i=0;i<ni;++i)
      for(int j=0;j<nj;++j)
        U[i][j] = PrimitiveToConserved(V[i][j]);

    ApplyLimitsToConserved("after V→U conversion", U);
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

/// Compute face lengths and outward normals on a structured curvilinear mesh.
/// - X, Y: node coordinates sized (imax+1)×(jmax+1)
/// - A_face_i, nx_face_i, ny_face_i: sized (imax+1)×jmax for vertical faces
/// - A_face_j, nx_face_j, ny_face_j: sized imax×(jmax+1) for horizontal faces
void computeMeshGeometry(
    const std::vector<std::vector<double>>& X,
    const std::vector<std::vector<double>>& Y,
    std::vector<std::vector<double>>& A_face_i,
    std::vector<std::vector<double>>& A_face_j,
    std::vector<std::vector<double>>& nx_face_i,
    std::vector<std::vector<double>>& ny_face_i,
    std::vector<std::vector<double>>& nx_face_j,
    std::vector<std::vector<double>>& ny_face_j
) {
    // Vertical (i‐faces): between node (i,j) → (i,j+1)
    for (int j = 0; j < jmax; ++j) {
        for (int i = 0; i <= imax; ++i) {
            double dx = X[i][j+1] - X[i][j];
            double dy = Y[i][j+1] - Y[i][j];
            double L  = std::hypot(dx, dy);
            A_face_i[i][j] = L;
            // outward normal (to the left cell): rotate (dx,dy) CCW 90° → ( dy, -dx )
            nx_face_i[i][j] =  dy / L;
            ny_face_i[i][j] = -dx / L;
        }
    }

    // Horizontal (j‐faces): between node (i,j) → (i+1,j)
    for (int j = 0; j <= jmax; ++j) {
        for (int i = 0; i < imax; ++i) {
            double dx = X[i+1][j] - X[i][j];
            double dy = Y[i+1][j] - Y[i][j];
            double L  = std::hypot(dx, dy);
            A_face_j[i][j] = L;
            // outward normal (to the bottom cell): rotate (dx,dy) CW 90° → ( -dy, dx )
            nx_face_j[i][j] = -dy / L;
            ny_face_j[i][j] =  dx / L;
        }
    }
}


// mirror-pad an im×jm array into (im+2*ghost)×(jm+2*ghost)
template<typename T>
static std::vector<std::vector<T>> padWithGhosts2D(
    const std::vector<std::vector<T>>& in,
    int ghost
) {
    int im = in.size();
    int jm = in.empty() ? 0 : in[0].size();
    int I  = im + 2*ghost;
    int J  = jm + 2*ghost;
    std::vector<std::vector<T>> out(I, std::vector<T>(J));
    // copy interior
    for(int j=0;j<jm;++j)
      for(int i=0;i<im;++i)
        out[i+ghost][j+ghost] = in[i][j];
    // mirror left/right
    for(int j=ghost;j<ghost+jm;++j)
      for(int g=0;g<ghost;++g){
        out[g][j]           = out[2*ghost-1-g][j];
        out[ghost+im+g][j]  = out[ghost+im-1-g][j];
      }
    // mirror bottom/top
    for(int i=0;i<I;++i)
      for(int g=0;g<ghost;++g){
        out[i][g]           = out[i][2*ghost-1-g];
        out[i][ghost+jm+g]  = out[i][ghost+jm-1-g];
      }
    return out;
}

/// Reads your .grd, builds Xn/Yn→cell-centers & face normals, pads ghosts.
/// Fills globals x_node,y_node and the passed-in eight arrays + returns bbox & dx,dy.
void readCurviMeshFromFile(
    const std::string & filepath,
    std::vector<std::vector<double>>& x_cell,
    std::vector<std::vector<double>>& y_cell,
    std::vector<std::vector<double>>& A_face_i,
    std::vector<std::vector<double>>& A_face_j,
    std::vector<std::vector<double>>& nx_face_i,
    std::vector<std::vector<double>>& ny_face_i,
    std::vector<std::vector<double>>& nx_face_j,
    std::vector<std::vector<double>>& ny_face_j,
    double & xmin, double & xmax,
    double & ymin, double & ymax,
    double & dx,   double & dy
) {
    std::ifstream in(filepath);
    if (!in) { std::cerr<<"[ERROR] Open "<<filepath<<"\n"; std::exit(1); }
    int nz; in>>nz>>imax>>jmax>>kmax;
    assert(nz==1 && kmax==2);

    // 1) Read node planes (keep k=0 only)
    std::vector<std::vector<double>> Xn(imax+1, std::vector<double>(jmax+1));
    std::vector<std::vector<double>> Yn(imax+1, std::vector<double>(jmax+1));
    double tmp;
    for(int k=0;k<kmax;++k)
      for(int j=0;j<=jmax;++j)
        for(int i=0;i<=imax;++i){
          in>>tmp;
          if(k==0) Xn[i][j]=tmp;
        }
    for(int k=0;k<kmax;++k)
      for(int j=0;j<=jmax;++j)
        for(int i=0;i<=imax;++i){
          in>>tmp;
          if(k==0) Yn[i][j]=tmp;
        }
    // skip Z
    for(int k=0;k<kmax;++k)
      for(int j=0;j<=jmax;++j)
        for(int i=0;i<=imax;++i) in>>tmp;

    // save for Tecplot
    x_node = Xn;
    y_node = Yn;

    // 2) bbox & estimates
    xmin = Xn[0][0];        xmax = Xn[imax][jmax];
    ymin = Yn[0][0];        ymax = Yn[imax][jmax];
    dx   = (xmax - xmin)/imax;
    dy   = (ymax - ymin)/jmax;

    // 3) cell‐centers
    x_cell.assign(imax, std::vector<double>(jmax));
    y_cell.assign(imax, std::vector<double>(jmax));
    for(int j=0;j<jmax;++j)
      for(int i=0;i<imax;++i){
        x_cell[i][j] = 0.25*(Xn[i][j]   + Xn[i+1][j]
                           + Xn[i][j+1] + Xn[i+1][j+1]);
        y_cell[i][j] = 0.25*(Yn[i][j]   + Yn[i+1][j]
                           + Yn[i][j+1] + Yn[i+1][j+1]);
      }

    // 4) face normals & lengths
    A_face_i .assign(imax+1, std::vector<double>(jmax));
    nx_face_i.assign(imax+1, std::vector<double>(jmax));
    ny_face_i.assign(imax+1, std::vector<double>(jmax));
    A_face_j .assign(imax,   std::vector<double>(jmax+1));
    nx_face_j.assign(imax,   std::vector<double>(jmax+1));
    ny_face_j.assign(imax,   std::vector<double>(jmax+1));
    computeMeshGeometry(
      Xn, Yn,
      A_face_i, A_face_j,
      nx_face_i, ny_face_i,
      nx_face_j, ny_face_j
    );

    // 5) now pad ghosts
    x_cell     = padWithGhosts2D(x_cell,     ghost);
    y_cell     = padWithGhosts2D(y_cell,     ghost);
    A_face_i   = padWithGhosts2D(A_face_i,   ghost);
    nx_face_i  = padWithGhosts2D(nx_face_i,  ghost);
    ny_face_i  = padWithGhosts2D(ny_face_i,  ghost);
    A_face_j   = padWithGhosts2D(A_face_j,   ghost);
    nx_face_j  = padWithGhosts2D(nx_face_j,  ghost);
    ny_face_j  = padWithGhosts2D(ny_face_j,  ghost);
}



// 2) Overwrite everything (cell-centers, face lengths, normals) with a uniform
//    Cartesian grid of size imax x jmax over [0,L]x[0,L].  This is debug mode.
//    After calling this, you can run your solver exactly as if you had a
//    curvilinear grid-because the same eight arrays are filled, but in a trivial
//    Cartesian way.
inline std::tuple<double,double,double,double,double,double>
convertToCartesianDebug(
    double L,
    std::vector<std::vector<double>>& x_cell,
    std::vector<std::vector<double>>& y_cell,
    std::vector<std::vector<double>>& A_face_i,
    std::vector<std::vector<double>>& A_face_j,
    std::vector<std::vector<double>>& nx_face_i,
    std::vector<std::vector<double>>& ny_face_i,
    std::vector<std::vector<double>>& nx_face_j,
    std::vector<std::vector<double>>& ny_face_j
) {
    int Ni = imax + 2 * ghost;
    int Nj = jmax + 2 * ghost;

    double dx = L / double(imax);  // imax interior cells => imax divisions
    double dy = L / double(jmax);

    // Allocate
    x_cell.assign(Ni, std::vector<double>(Nj));
    y_cell.assign(Ni, std::vector<double>(Nj));

    for (int i = 0; i < Ni; ++i) {
        for (int j = 0; j < Nj; ++j) {
            x_cell[i][j] = dx * (i - ghost + 0.5);  // center of each cell
            y_cell[i][j] = dy * (j - ghost + 0.5);
        }
    }

    // Face areas and normals
    A_face_i.assign(Ni + 1, std::vector<double>(Nj, dy));
    nx_face_i.assign(Ni + 1, std::vector<double>(Nj, 1.0));
    ny_face_i.assign(Ni + 1, std::vector<double>(Nj, 0.0));

    A_face_j.assign(Ni, std::vector<double>(Nj + 1, dx));
    nx_face_j.assign(Ni, std::vector<double>(Nj + 1, 0.0));
    ny_face_j.assign(Ni, std::vector<double>(Nj + 1, 1.0));

    // Physical domain
    double xmin = 0.0;
    double xmax = L;
    double ymin = 0.0;
    double ymax = L;

    return {xmin, xmax, ymin, ymax, dx, dy};
}




//----------------------------------------------------------------------
// (A) Define supersonic and subsonic MMS constants
//----------------------------------------------------------------------

// Constants for MMS field construction
struct MmsParams {
    double rho0,   rho_x,   rho_y, rho_z, a_rho_x, a_rho_y, a_rho_z;
    double u0,     u_x,     u_y, u_z, a_u_x,   a_u_y, a_u_z;
    double v0,     v_x,     v_y, v_z, a_v_x,   a_v_y, a_v_z;
    double p0,     p_x,     p_y, p_z, a_p_x,   a_p_y, a_p_z;
};

// Supersonic constants (MMS case 1)
constexpr MmsParams mmsSup = {
    1.0,   0.15,  -0.10, 0,  1.0,    0.50, 0,       // rho
    800.0, 50.0, -30.0, 0,  1.5,    0.60, 0,       // u
    800.0, -75.0, 40.0, 0,   0.5,    (2.0 / 3.0), 0,// v
    100000.0, 20000.0, 50000.0, 0, 2.0, 1.0, 0     // p
};

// Subsonic constants (MMS case 2)
constexpr MmsParams mmsSub = {
    1.0,   0.15,  -0.10, 0,  1.0,   0.50, 0,       // rho
    70.0,  5.0,   -7.0, 0,  1.5,   0.60, 0,      // u
    90.0, -15.0,   8.5, 0,  0.5,   (2.0 / 3.0), 0,// v
    100000.0, 20000.0, 50000.0, 0, 2.0, 1.0, 0     // p
};

constexpr double PI = 3.14159265358979323846;


inline double rho_mms(int mmsCase, double L, double x, double y) {
    assert(mmsCase == 1 || mmsCase == 2);
    const auto& C = (mmsCase == 1 ? mmsSup : mmsSub);

    // Fortran: rho0 + rhoy*cos((pi*y)/(2*length)) + rhox*sin((pi*x)/length)
    return C.rho0
         + C.rho_y * std::cos((PI * y) / (2.0 * L))
         + C.rho_x * std::sin((PI * x) / L);
}

inline double uvel_mms(int mmsCase, double L, double x, double y) {
    assert(mmsCase == 1 || mmsCase == 2);
    const auto& C = (mmsCase == 1 ? mmsSup : mmsSub);

    return C.u0
         + C.u_y * std::cos((3.0 * PI * y) / (5.0 * L))
         + C.u_x * std::sin((3.0 * PI * x) / (2.0 * L));
}

inline double vvel_mms(int mmsCase, double L, double x, double y) {
    assert(mmsCase == 1 || mmsCase == 2);
    const auto& C = (mmsCase == 1 ? mmsSup : mmsSub);

    // Fortran:
    //   vvel0
    // + vvelx * cos((pi*x)/(two*length))
    // + vvely * sin((two*pi*y)/(three*length))
    return C.v0
         + C.v_x * std::cos((PI * x) / (2.0 * L))
         + C.v_y * std::sin((2.0 * PI * y) / (3.0 * L));
}

inline double press_mms(int mmsCase, double L, double x, double y) {
    assert(mmsCase == 1 || mmsCase == 2);
    const auto& C = (mmsCase == 1 ? mmsSup : mmsSub);

    // Fortran: press0
    //       + pressx * cos((two*pi*x)/length)
    //       + pressy * sin((pi*y)/length)
    return C.p0
         + C.p_x * std::cos((2.0 * PI * x) / L)
         + C.p_y * std::sin((PI * y) / L);
}

//----------------------------------------------------------------------
// (C) Exact “source” fields (mass‐equation, x‐momentum, y‐momentum, energy)
//   Each formula is copied verbatim from your Fortran → C++ conversion,
//   but replacing every “rho0, rhox, …” with C.<name> based on mmsCase.
//----------------------------------------------------------------------

inline double rmassconv(int mmsCase,
                        double L,
                        double x,
                        double y)
{
    assert(mmsCase==1 || mmsCase==2);
    const auto& C = (mmsCase==1 ? mmsSup : mmsSub);

    return 
      // term1
      (3.0*PI * C.u_x 
       * std::cos((3.0*PI*x)/(2.0*L))
       * (C.rho0 
          + C.rho_y*std::cos((PI*y)/(2.0*L)) 
          + C.rho_x*std::sin((PI*x)/(L))))
      /(2.0*L)

      // + term2
      + (2.0*PI * C.v_y 
         * std::cos((2.0*PI*y)/(3.0*L))
         * (C.rho0 
            + C.rho_y*std::cos((PI*y)/(2.0*L)) 
            + C.rho_x*std::sin((PI*x)/(L))))
      /(3.0*L)

      // + term3
      + (PI * C.rho_x 
         * std::cos((PI*x)/(L))
         * (C.u0 
            + C.u_y*std::cos((3.0*PI*y)/(5.0*L)) 
            + C.u_x*std::sin((3.0*PI*x)/(2.0*L))))
      /(L)

      // – term4
      - (PI * C.rho_y 
         * std::sin((PI*y)/(2.0*L))
         * (C.v0 
            + C.v_x*std::cos((PI*x)/(2.0*L)) 
            + C.v_y*std::sin((2.0*PI*y)/(3.0*L))))
      /(2.0*L);
}

inline double xmtmconv(int mmsCase, double L, double x, double y) {
    assert(mmsCase==1 || mmsCase==2);
    const auto& C = (mmsCase==1 ? mmsSup : mmsSub);

    // term1
    double term1 = (3.0*PI * C.u_x
                   * std::cos((3.0*PI*x)/(2.0*L))
                   * (C.rho0
                      + C.rho_y * std::cos((PI*y)/(2.0*L))
                      + C.rho_x * std::sin((PI*x)/L))
                   * (C.u0
                      + C.u_y * std::cos((3.0*PI*y)/(5.0*L))
                      + C.u_x * std::sin((3.0*PI*x)/(2.0*L))))
                  / L;

    // term2
    double term2 = (2.0*PI * C.v_y
                   * std::cos((2.0*PI*y)/(3.0*L))
                   * (C.rho0
                      + C.rho_y * std::cos((PI*y)/(2.0*L))
                      + C.rho_x * std::sin((PI*x)/L))
                   * (C.u0
                      + C.u_y * std::cos((3.0*PI*y)/(5.0*L))
                      + C.u_x * std::sin((3.0*PI*x)/(2.0*L))))
                  / (3.0*L);

    // u‐piece reused in term3 & term5
    double uPiece = (C.u0
                   + C.u_y * std::cos((3.0*PI*y)/(5.0*L))
                   + C.u_x * std::sin((3.0*PI*x)/(2.0*L)));

    // term3
    double term3 = (PI * C.rho_x
                   * std::cos((PI*x)/L)
                   * (uPiece * uPiece))
                  / L;

    // term4
    double term4 = (2.0*PI * C.p_x * std::sin((2.0*PI*x)/L))
                   / L;

    // v‐piece reused in term5 & term6
    double vPiece = (C.v0
                   + C.v_x * std::cos((PI*x)/(2.0*L))
                   + C.v_y * std::sin((2.0*PI*y)/(3.0*L)));

    // term5
    double term5 = (PI * C.rho_y
                   * uPiece
                   * std::sin((PI*y)/(2.0*L))
                   * vPiece)
                  / (2.0*L);

    // rho‐piece reused in term6
    double rhoPiece = (C.rho0
                     + C.rho_y * std::cos((PI*y)/(2.0*L))
                     + C.rho_x * std::sin((PI*x)/L));

    // term6
    double term6 = (3.0*PI * C.u_y
                   * rhoPiece
                   * std::sin((3.0*PI*y)/(5.0*L))
                   * vPiece)
                  / (5.0*L);

    return term1
         + term2
         + term3
         - term4
         - term5
         - term6;
}

inline double ymtmconv(int mmsCase, double L, double x, double y) {
    assert(mmsCase==1 || mmsCase==2);
    const auto& C = (mmsCase==1 ? mmsSup : mmsSub);

    // reusable pieces
    double rhoPiece = (C.rho0
                     + C.rho_y * std::cos((PI * y)/(2.0*L))
                     + C.rho_x * std::sin((PI * x)/L));
    double uPiece   = (C.u0
                     + C.u_y * std::cos((3.0*PI * y)/(5.0*L))
                     + C.u_x * std::sin((3.0*PI * x)/(2.0*L)));
    double vPiece   = (C.v0
                     + C.v_x * std::cos((PI * x)/(2.0*L))
                     + C.v_y * std::sin((2.0*PI * y)/(3.0*L)));

    // term1:  (pi*pressy*cos((pi*y)/L)) / L
    double term1 = (PI * C.p_y * std::cos((PI * y)/L)) / L;

    // term2:  (pi*uvelx*sin((pi*x)/(2L)) * rhoPiece * uPiece) / (2L)
    double term2 = (PI * C.u_x
                  * std::sin((PI * x)/(2.0*L))
                  * rhoPiece * uPiece)
                  / (2.0*L);

    // term3:  (3*pi*uvelx*cos((3*pi*x)/(2L)) * rhoPiece * vPiece) / (2L)
    double term3 = (3.0*PI * C.u_x
                  * std::cos((3.0*PI * x)/(2.0*L))
                  * rhoPiece * vPiece)
                  / (2.0*L);

    // term4:  (4*pi*vvely*cos((2*pi*y)/(3L)) * rhoPiece * vPiece) / (3L)
    double term4 = (4.0*PI * C.v_y
                  * std::cos((2.0*PI * y)/(3.0*L))
                  * rhoPiece * vPiece)
                  / (3.0*L);

    // term5:  (pi*rhox*cos((pi*x)/L) * uPiece * vPiece) / L
    double term5 = (PI * C.rho_x
                  * std::cos((PI * x)/L)
                  * uPiece * vPiece)
                  / L;

    // term6:  (pi*rhoy*sin((pi*y)/(2L)) * vPiece^2) / (2L)
    double term6 = (PI * C.rho_y
                  * std::sin((PI * y)/(2.0*L))
                  * vPiece * vPiece)
                  / (2.0*L);

    // assemble: +term1 -term2 +term3 +term4 +term5 -term6
    return term1
         - term2
         + term3
         + term4
         + term5
         - term6;
}

inline double energyconv(int mmsCase,
                         double gamma,
                         double L,
                         double x,
                         double y)
{
    assert(mmsCase == 1 || mmsCase == 2);
    const auto& C = (mmsCase == 1 ? mmsSup : mmsSub);

    // ρ field
    double rhoPhi = C.rho0
                  + C.rho_y * std::cos((PI * y) / (2.0 * L))
                  + C.rho_x * std::sin((PI * x) / L);

    // u‐velocity field
    double uPhi = C.u0
                + C.u_y * std::cos((3.0 * PI * y) / (5.0 * L))
                + C.u_x * std::sin((3.0 * PI * x) / (2.0 * L));

    // v‐velocity field
    double vPhi = C.v0
                + C.v_x * std::cos((PI * x) / (2.0 * L))
                + C.v_y * std::sin((2.0 * PI * y) / (3.0 * L));

    // pressure field
    double pPhi = C.p0
                + C.p_x * std::cos((2.0 * PI * x) / L)
                + C.p_y * std::sin((PI * y) / L);

    // specific total energy
    double vel2 = uPhi*uPhi + vPhi*vPhi;
    double Eeq  = 0.5*vel2 + pPhi/((gamma - 1.0)*rhoPhi);

    // bracket1
    double A = -2.0*PI * C.p_x
             * std::sin((2.0*PI * x) / L)
             / L;
    double B = rhoPhi * (
                   -2.0*PI * C.p_x
                    * std::sin((2.0*PI * x) / L)
                    / ((gamma - 1.0)*L*rhoPhi)
                   + ((3.0*PI * C.u_x
                       * std::cos((3.0*PI * x)/(2.0*L))
                       * uPhi
                     - PI    * C.v_x
                       * std::sin((PI * x)/(2.0*L))
                       * vPhi)
                      / L)
               ) * 0.5;
    double C1 = -PI * C.rho_x
              * std::cos((PI * x) / L)
              * pPhi
              / ((gamma - 1.0)*L*rhoPhi*rhoPhi);
    double bracket1 = A + B + C1;

    // terms 1–6
    double term1 = uPhi * bracket1;

    double term2 = (PI * C.rho_x
                   * std::cos((PI*x)/L)
                   * Eeq)
                  / L;

    double term3 = (3.0*PI * C.u_x
                   * std::cos((3.0*PI*x)/(2.0*L))
                   * (pPhi + rhoPhi*Eeq))
                  / (2.0*L);

    double term4 = (2.0*PI * C.v_y
                   * std::cos((2.0*PI*y)/(3.0*L))
                   * (pPhi + rhoPhi*Eeq))
                  / (3.0*L);

    double bracket2 = (PI * C.p_y * std::cos((PI*y)/L)) / L
                    - (PI * C.rho_y * std::sin((PI*y)/(2.0*L)) * vPhi) / (2.0*L);
    double term5 = vPhi * bracket2;

    double bracket3 = (PI * C.p_y * std::cos((PI*y)/L)) / ((gamma - 1.0)*L*rhoPhi)
                    + ((-6.0*PI * C.u_y * uPhi * std::sin((3.0*PI*y)/(5.0*L)))
                       / (5.0*L)
                       + (4.0*PI * C.v_y * std::cos((2.0*PI*y)/(3.0*L)) * vPhi)
                         / (3.0*L)) * 0.5
                    + (PI * C.rho_y * std::sin((PI*y)/(2.0*L)) * pPhi)
                      / (2.0*(gamma - 1.0)*L*rhoPhi*rhoPhi);
    double term6 = rhoPhi * bracket3;

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
    int i , int j ,
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
    std::ofstream file(filename, std::ios::app);
    if (!file) {
        std::cerr << "Error: cannot open " << filename << " for writing\n";
        return;
    }

    file << "TITLE = \"2D Euler MMS Solution\"\n";
    file << "VARIABLES = \"X\" \"Y\" \"rho\" \"u\" \"v\" \"P\" \"Mach\"\n";
    file << "ZONE I=" << imax+1 << " J=" << jmax+1
         << ", DATAPACKING=BLOCK\n";
    file << "VARLOCATION=([3-7]=CELLCENTERED)\n";

    // Helper to guard against NaNs
    auto safe = [](double v) {
        return std::isnan(v) ? -999.9 : v;
    };

    // 3. Dump X block (NODAL)
    for (int j = 0; j <= jmax; ++j)
      for (int i = 0; i <= imax; ++i)
        file << x_node[i][j] << "\n";


    // 4. Dump Y block (NODAL)
    for (int j = 0; j <= jmax; ++j)
      for (int i = 0; i <= imax; ++i)
        file << y_node[i][j] << "\n";


    for (int j = 0; j < jmax; ++j)
      for (int i = 0; i < imax; ++i)
        file << V[i+ghost][j+ghost].rho << "\n";

    for (int j = 0; j < jmax; ++j)
      for (int i = 0; i < imax; ++i)
        file << V[i+ghost][j+ghost].u << "\n";

    for (int j = 0; j < jmax; ++j)
      for (int i = 0; i < imax; ++i)
        file << V[i+ghost][j+ghost].v << "\n";

    for (int j = 0; j < jmax; ++j)
      for (int i = 0; i < imax; ++i)
        file << V[i+ghost][j+ghost].P << "\n";

    for (int j = 0; j < jmax; ++j)
      for (int i = 0; i < imax; ++i) {
        double rho = V[i+ghost][j+ghost].rho;
        double u = V[i+ghost][j+ghost].u;
        double v = V[i+ghost][j+ghost].v;
        double P = V[i+ghost][j+ghost].P;
        double a = std::sqrt(gamma * P / rho);
        double mach = std::sqrt(u*u + v*v) / a;
        file << mach << "\n";
    }
    file.close();
    std::cout << "[INFO] Appended 2D solution for iter=" << iter
              << " to " << filename << " in BLOCK format\n";
}



int main() {

  int fluxOrder = 2;  // Second-order MUSCL
  double kappa = 0.5; // Choose limiter parameter (e.g., -1, 0, 0.5, or 1)
  bool freezeLimiter = true;  // Set to true for pure MUSCL (no limiter)

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
    const std::string meshFile = R"(C:\Users\monicashanmugam\OneDrive - Virginia Tech\Desktop\Virginia Tech\CFD\Projectv2\Project_Files\Project_Files\Grids\curviliniear-grids\curv2d17.grd)";

    // 2) Decide whether you want to run in DEBUG (Cartesian) mode
    bool debugMode = false;  // set true if you want a simple Cartesian mesh  
    
    double xmin, xmax, ymin, ymax, dx, dy;
    
    if (!debugMode) {
      readCurviMeshFromFile(meshFile,x_cell, y_cell,A_face_i, A_face_j,nx_face_i, ny_face_i,nx_face_j, ny_face_j,xmin, xmax, ymin, ymax,dx, dy);

      std::cout << "[INFO] Loaded curvi mesh with imax = " << imax << ", jmax = " << jmax << "\n";
    }
    else {
      // in debug mode, only read imax/jmax and then overwrite with Cartesian:
      
        std::ifstream in(meshFile);
        if(!in){
          std::cerr << "Error cannot open mesh file. \n";
          std::exit(EXIT_FAILURE);
        }
        int nz, kmax;
        in >> nz >> imax >> jmax >> kmax;
        assert(nz == 1 && kmax == 2);
    
      std::tie(xmin, xmax, ymin, ymax, dx, dy) = convertToCartesianDebug(L, x_cell, y_cell, A_face_i, A_face_j, nx_face_i, ny_face_i, nx_face_j, ny_face_j);
      std::cout << "[INFO] Using debug Cartesian mesh: imax="<<imax<<", jmax="<<jmax<<"\n";

      // 🟩 Add this to generate nodal coordinates
      x_node.resize(imax + 1, std::vector<double>(jmax + 1));
      y_node.resize(imax + 1, std::vector<double>(jmax + 1));

      for (int j = 0; j <= jmax; ++j) {
        for (int i = 0; i <= imax; ++i) {
            x_node[i][j] = i * dx;
            y_node[i][j] = j * dy;
        }
      
      }

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

    //     // right after initializeMMS + BCs:
    // OutputSolution2D(solFile, 0);


    std::vector<std::vector<Conserved>> S;
    computeSourceTermsMMS(mmsCase, gamma, L, x_cell, y_cell, S);

    GlobalPrimitiveToConserved();

    double dt = computeTimeStep(V, cellVolume, A_face_i, A_face_j, nx_face_i, ny_face_i, nx_face_j, ny_face_j);
    std::cout << "[INFO] Computed time step: dt = " << dt << "\n";

    

    std::vector<std::vector<Conserved>> R_int(imax, std::vector<Conserved>(jmax));  // Residuals

        // Pseudo-time marching parameters
    const int maxIter = 1;
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
