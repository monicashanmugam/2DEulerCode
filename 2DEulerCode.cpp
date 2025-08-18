#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <filesystem>  
#include <tuple>
#include <cstdlib>



using namespace std;

namespace fs = std::filesystem;

//-----------------------------------------------------
// Global constants and parameters
//-----------------------------------------------------
constexpr double Rgas      = 287.0;  // [J/(kg·K)]
constexpr double gamma = 1.4;
constexpr double L = 1.0; //domain length
static constexpr int ghost = 2;
const double CFL = 0.02;    // CFL number for time step control
const double epsM = 0.01;  // Minimum Mach number allowed
const double tolerance = 1e-12;
const double delta = 1e-6; // For computing the r

std::vector<std::vector<double>> x_node;
std::vector<std::vector<double>> y_node;

static_assert(ghost >= 2, "MUSCL reconstruction needs at least 2 ghost cells.");



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

// ---------- tiny helpers ----------
inline double tri_area(double xa,double ya,double xb,double yb,double xc,double yc){
    return 0.5*std::abs((xb-xa)*(yc-ya) - (yb-ya)*(xc-xa));
}

template<typename T>
std::vector<std::vector<T>> padWithGhosts2D(const std::vector<std::vector<T>>& in, int g){
    int im = (int)in.size(); if (!im) return in;
    int jm = (int)in[0].size();
    std::vector<std::vector<T>> out(im+2*g, std::vector<T>(jm+2*g));
    for(int i=0;i<im;++i) for(int j=0;j<jm;++j) out[i+g][j+g]=in[i][j];
    for(int j=g;j<g+jm;++j) for(int q=0;q<g;++q){ out[g-1-q][j]=out[g+q][j]; out[g+im+q][j]=out[g+im-1-q][j]; }
    for(int i=0;i<im+2*g;++i) for(int q=0;q<g;++q){ out[i][g-1-q]=out[i][g+q]; out[i][g+jm+q]=out[i][g+jm-1-q]; }
    return out;
}
static void padFacesScalar(std::vector<std::vector<double>>& A, int g){
    int I=(int)A.size(), J=(I? (int)A[0].size():0); if(!I||!J||!g) return;
    for(int j=g;j<=J-g-1;++j) for(int q=0;q<g;++q){ A[q][j]=A[2*g-1-q][j]; A[I-g+q][j]=A[I-g-1-q][j]; }
    for(int i=0;i<I;++i) for(int q=0;q<g;++q){ A[i][q]=A[i][2*g-1-q]; A[i][J-g+q]=A[i][J-g-1-q]; }
}
static void padFacesNormal(std::vector<std::vector<double>>& nx,
                           std::vector<std::vector<double>>& ny, int g){
    int I=(int)nx.size(), J=(I? (int)nx[0].size():0); if(!I||!J||!g) return;
    for(int j=g;j<=J-g-1;++j) for(int q=0;q<g;++q){
        int il=2*g-1-q, ir=I-g-1-q;
        nx[q][j]=-nx[il][j]; ny[q][j]=-ny[il][j];
        nx[I-g+q][j]=-nx[ir][j]; ny[I-g+q][j]=-ny[ir][j];
    }
    for(int i=0;i<I;++i) for(int q=0;q<g;++q){
        int jb=2*g-1-q, jt=J-g-1-q;
        nx[i][q]=-nx[i][jb]; ny[i][q]=-ny[i][jb];
        nx[i][J-g+q]=-nx[i][jt]; ny[i][J-g+q]=-ny[i][jt];
    }
}
// unsnake rows/cols so adjacency is consistent
static void harmonize_ordering(std::vector<std::vector<double>>& Xn,
                               std::vector<std::vector<double>>& Yn){
    int Ni=(int)Xn.size(), Nj=(int)Xn[0].size();
    auto rowsign=[&](int j){ return (Xn[Ni-1][j]-Xn[0][j] >= 0) ? 1: -1; };
    auto colsign=[&](int i){ return (Yn[i][Nj-1]-Yn[i][0] >= 0) ? 1: -1; };
    int rs=rowsign(0); for(int j=1;j<Nj;++j) if (rowsign(j)*rs<0)
        for(int i=0;i<Ni/2;++i){ std::swap(Xn[i][j],Xn[Ni-1-i][j]); std::swap(Yn[i][j],Yn[Ni-1-i][j]); }
    int cs=colsign(0); for(int i=1;i<Ni;++i) if (colsign(i)*cs<0)
        for(int j=0;j<Nj/2;++j){ std::swap(Xn[i][j],Xn[i][Nj-1-j]); std::swap(Yn[i][j],Yn[i][Nj-1-j]); }
}

// Reads a 2D curvilinear .grd with header:
//   nz
//   Ni_nodes  Nj_nodes  kplanes
// X-block (kplanes planes; use plane 0), Y-block (kplanes; use plane 0), Z-block (ignored)
// Builds: x_node/y_node (no ghosts), and padded x_cell/y_cell, faces, normals, cellVolume.
// Sets global imax = Ni-1, jmax = Nj-1. Uses global `ghost`.
void readCurviMeshFromFile(
    const std::string& path,
    std::vector<std::vector<double>>& x_cell,
    std::vector<std::vector<double>>& y_cell,
    std::vector<std::vector<double>>& A_face_i,
    std::vector<std::vector<double>>& A_face_j,
    std::vector<std::vector<double>>& nx_face_i,
    std::vector<std::vector<double>>& ny_face_i,
    std::vector<std::vector<double>>& nx_face_j,
    std::vector<std::vector<double>>& ny_face_j,
    std::vector<std::vector<double>>& cellVolume, // <-- padded volume
    double& xmin, double& xmax, double& ymin, double& ymax, double& dx, double& dy
) {
    std::ifstream in(path);
    if (!in) { std::cerr << "[ERROR] Cannot open mesh file: " << path << "\n"; std::exit(EXIT_FAILURE); }

    int nz, Ni_nodes, Nj_nodes, kplanes;
    in >> nz >> Ni_nodes >> Nj_nodes >> kplanes;
    if (!(nz == 1 && kplanes >= 1)) {
        std::cerr << "[ERROR] Bad .grd header (nz or kplanes)\n"; std::exit(EXIT_FAILURE);
    }

    // cells = nodes - 1
    imax = Ni_nodes - 1;
    jmax = Nj_nodes - 1;

    // --- read X (take plane 0) ---
    std::vector<std::vector<double>> Xn(Ni_nodes, std::vector<double>(Nj_nodes));
    std::vector<std::vector<double>> Yn(Ni_nodes, std::vector<double>(Nj_nodes));
    double tmp;
    for (int k = 0; k < kplanes; ++k)
        for (int j = 0; j < Nj_nodes; ++j)
            for (int i = 0; i < Ni_nodes; ++i) { in >> tmp; if (k == 0) Xn[i][j] = tmp; }
    // --- read Y (take plane 0) ---
    for (int k = 0; k < kplanes; ++k)
        for (int j = 0; j < Nj_nodes; ++j)
            for (int i = 0; i < Ni_nodes; ++i) { in >> tmp; if (k == 0) Yn[i][j] = tmp; }
    // --- skip Z ---
    for (int k = 0; k < kplanes; ++k)
        for (int j = 0; j < Nj_nodes; ++j)
            for (int i = 0; i < Ni_nodes; ++i) in >> tmp;

    // --- harmonize ordering (unsnake rows/cols) ---
    auto reverse_row = [&](int j) {
        for (int i = 0; i < Ni_nodes/2; ++i) { std::swap(Xn[i][j], Xn[Ni_nodes-1-i][j]); std::swap(Yn[i][j], Yn[Ni_nodes-1-i][j]); }
    };
    auto reverse_col = [&](int i) {
        for (int j = 0; j < Nj_nodes/2; ++j) { std::swap(Xn[i][j], Xn[i][Nj_nodes-1-j]); std::swap(Yn[i][j], Yn[i][Nj_nodes-1-j]); }
    };
    auto rowsign = [&](int j){ return (Xn[Ni_nodes-1][j] - Xn[0][j]) >= 0 ? 1 : -1; };
    auto colsign = [&](int i){ return (Yn[i][Nj_nodes-1] - Yn[i][0]) >= 0 ? 1 : -1; };
    int rs = rowsign(0); for (int j = 1; j < Nj_nodes; ++j) if (rowsign(j)*rs < 0) reverse_row(j);
    int cs = colsign(0); for (int i = 1; i < Ni_nodes; ++i) if (colsign(i)*cs < 0) reverse_col(i);

    // --- save nodal arrays for Tecplot (no ghosts) ---
    x_node = Xn;  y_node = Yn;

    // --- centers (interior imax x jmax) ---
    x_cell.assign(imax, std::vector<double>(jmax));
    y_cell.assign(imax, std::vector<double>(jmax));
    for (int j = 0; j < jmax; ++j)
        for (int i = 0; i < imax; ++i) {
            x_cell[i][j] = 0.25*(Xn[i][j] + Xn[i+1][j] + Xn[i][j+1] + Xn[i+1][j+1]);
            y_cell[i][j] = 0.25*(Yn[i][j] + Yn[i+1][j] + Yn[i][j+1] + Yn[i+1][j+1]);
        }

    // --- face lengths & outward normals (interior faces) ---
    A_face_i.assign(imax+1, std::vector<double>(jmax));
    nx_face_i.assign(imax+1, std::vector<double>(jmax));
    ny_face_i.assign(imax+1, std::vector<double>(jmax));
    for (int j = 0; j < jmax; ++j)
        for (int i = 0; i <= imax; ++i) {
            double dx_ = Xn[i][j+1] - Xn[i][j];
            double dy_ = Yn[i][j+1] - Yn[i][j];
            double L = std::hypot(dx_, dy_);
            A_face_i[i][j] = L;
            nx_face_i[i][j] =  dy_ / L;  // outward from left cell
            ny_face_i[i][j] = -dx_ / L;
        }

    A_face_j.assign(imax, std::vector<double>(jmax+1));
    nx_face_j.assign(imax, std::vector<double>(jmax+1));
    ny_face_j.assign(imax, std::vector<double>(jmax+1));
    for (int j = 0; j <= jmax; ++j)
        for (int i = 0; i < imax; ++i) {
            double dx_ = Xn[i+1][j] - Xn[i][j];
            double dy_ = Yn[i+1][j] - Yn[i][j];
            double L = std::hypot(dx_, dy_);
            A_face_j[i][j] = L;
            nx_face_j[i][j] = -dy_ / L;  // outward from bottom cell
            ny_face_j[i][j] =  dx_ / L;
        }

    // --- cell areas from nodes (interior) ---
    cellVolume.assign(imax, std::vector<double>(jmax));
    auto tri_area = [](double xa,double ya,double xb,double yb,double xc,double yc){
        return 0.5*std::abs((xb-xa)*(yc-ya) - (yb-ya)*(xc-xa));
    };
    for (int j = 0; j < jmax; ++j)
        for (int i = 0; i < imax; ++i) {
            double x00 = Xn[i][j],     y00 = Yn[i][j];
            double x10 = Xn[i+1][j],   y10 = Yn[i+1][j];
            double x11 = Xn[i+1][j+1], y11 = Yn[i+1][j+1];
            double x01 = Xn[i][j+1],   y01 = Yn[i][j+1];
            cellVolume[i][j] = tri_area(x00,y00,x10,y10,x11,y11) + tri_area(x00,y00,x11,y11,x01,y01);
        }

    // --- padding helpers (mirror; normals flip sign in belts) ---
    auto padScalar = [&](const std::vector<std::vector<double>>& in)->std::vector<std::vector<double>>{
        int g = ghost;
        if (in.empty() || in[0].empty() || g<=0) return in;
        int I = (int)in.size(), J = (int)in[0].size();
        std::vector<std::vector<double>> out(I+2*g, std::vector<double>(J+2*g));
        for(int i=0;i<I;++i) for(int j=0;j<J;++j) out[i+g][j+g]=in[i][j];
        for(int j=g;j<g+J;++j) for(int q=0;q<g;++q){ out[g-1-q][j]=out[g+q][j]; out[g+I+q][j]=out[g+I-1-q][j]; }
        for(int i=0;i<I+2*g;++i) for(int q=0;q<g;++q){ out[i][g-1-q]=out[i][g+q]; out[i][g+J+q]=out[i][g+J-1-q]; }
        return out;
    };
    auto padNormals = [&](std::vector<std::vector<double>>& nx,
                          std::vector<std::vector<double>>& ny){
        int g = ghost; if (nx.empty()||nx[0].empty()||g<=0) return;
        int I = (int)nx.size(), J = (int)nx[0].size();
        // left/right belts
        for(int j=g;j<=J-g-1;++j) for(int q=0;q<g;++q){
            int il = 2*g-1-q, ir = I-g-1-q;
            nx[q][j]        = -nx[il][j]; ny[q][j]        = -ny[il][j];
            nx[I-g+q][j]    = -nx[ir][j]; ny[I-g+q][j]    = -ny[ir][j];
        }
        // bottom/top belts
        for(int i=0;i<I;++i) for(int q=0;q<g;++q){
            int jb = 2*g-1-q, jt = J-g-1-q;
            nx[i][q]        = -nx[i][jb]; ny[i][q]        = -ny[i][jb];
            nx[i][J-g+q]    = -nx[i][jt]; ny[i][J-g+q]    = -ny[i][jt];
        }
    };

    // --- pad everything used by the solver ---
    x_cell     = padScalar(x_cell);
    y_cell     = padScalar(y_cell);
    A_face_i   = padScalar(A_face_i);
    nx_face_i  = padScalar(nx_face_i);
    ny_face_i  = padScalar(ny_face_i);
    A_face_j   = padScalar(A_face_j);
    nx_face_j  = padScalar(nx_face_j);
    ny_face_j  = padScalar(ny_face_j);
    cellVolume = padScalar(cellVolume);

    // fix normals in ghost belts
    padNormals(nx_face_i, ny_face_i);
    padNormals(nx_face_j, ny_face_j);

    std::cout << "[READ] nodes: " << Ni_nodes << " x " << Nj_nodes
              << "  => cells: " << imax << " x " << jmax
              << "  (ghost=" << ghost << ")\n";
}




// Build a uniform Cartesian mesh over [0,Lx]×[0,Ly] with imax×jmax CELLS.
// Fills x_node/y_node (no ghosts) and padded x_cell/y_cell, faces, normals, cellVolume.
// Shapes (with ghost=g):
//  - V,U,x_cell,y_cell,cellVolume: (imax+2g) × (jmax+2g)
//  - A_face_i,nx_face_i,ny_face_i: (imax+1+2g) × (jmax+2g)
//  - A_face_j,nx_face_j,ny_face_j: (imax+2g)   × (jmax+1+2g)
void buildCartesianDebug(
    int imax, int jmax, double Lx, double Ly, int g,
    std::vector<std::vector<double>>& x_cell,
    std::vector<std::vector<double>>& y_cell,
    std::vector<std::vector<double>>& A_face_i,
    std::vector<std::vector<double>>& A_face_j,
    std::vector<std::vector<double>>& nx_face_i,
    std::vector<std::vector<double>>& ny_face_i,
    std::vector<std::vector<double>>& nx_face_j,
    std::vector<std::vector<double>>& ny_face_j,
    std::vector<std::vector<double>>& cellVolume,
    std::vector<std::vector<double>>& x_node,
    std::vector<std::vector<double>>& y_node,
    double &dx, double &dy
) {
    // Nodes (no ghosts)
    x_node.assign(imax+1, std::vector<double>(jmax+1));
    y_node.assign(imax+1, std::vector<double>(jmax+1));

    dx = Lx / double(imax);
    dy = Ly / double(jmax);

    for (int j=0; j<=jmax; ++j)
      for (int i=0; i<=imax; ++i) {
        x_node[i][j] = i*dx;
        y_node[i][j] = j*dy;
      }

    // Padded sizes
    const int Ni = imax + 2*g;
    const int Nj = jmax + 2*g;

    // Cell centers (padded)
    x_cell.assign(Ni, std::vector<double>(Nj));
    y_cell.assign(Ni, std::vector<double>(Nj));
    for (int j=0; j<Nj; ++j) {
      for (int i=0; i<Ni; ++i) {
        // physical interior is i∈[g..g+imax-1], j∈[g..g+jmax-1]
        const double xc = (i - g + 0.5)*dx;
        const double yc = (j - g + 0.5)*dy;
        x_cell[i][j] = xc;
        y_cell[i][j] = yc;
      }
    }

    // Faces (padded) — constant lengths, axis-aligned normals
    A_face_i.assign(Ni+1, std::vector<double>(Nj, dy));
    nx_face_i.assign(Ni+1, std::vector<double>(Nj, 1.0));
    ny_face_i.assign(Ni+1, std::vector<double>(Nj, 0.0));

    A_face_j.assign(Ni, std::vector<double>(Nj+1, dx));
    nx_face_j.assign(Ni, std::vector<double>(Nj+1, 0.0));
    ny_face_j.assign(Ni, std::vector<double>(Nj+1, 1.0));

    // Cell volumes (padded)
    cellVolume.assign(Ni, std::vector<double>(Nj, dx*dy));
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

inline double rmassconv(int mmsCase, double L, double x, double y) {
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

inline double energyconv(int mmsCase, double gamma, double L, double x, double y){
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

void initializeMMS(
    int mmsCase, double L,
    const std::vector<std::vector<double>>& x_cell,
    const std::vector<std::vector<double>>& y_cell
) {
    const int Ni = imax + 2*ghost;
    const int Nj = jmax + 2*ghost;

    // Ensure V is sized (if not already sized before this call)
    if ((int)V.size() != Ni || (int)V[0].size() != Nj) {
        V.assign(Ni, std::vector<Primitive>(Nj));
    }

    // Fill ONLY physical cells: i∈[ghost .. ghost+imax-1], j∈[ghost .. ghost+jmax-1]
    for (int i = ghost; i < ghost + imax; ++i) {
        for (int j = ghost; j < ghost + jmax; ++j) {
            double x = x_cell[i][j];
            double y = y_cell[i][j];

            V[i][j].rho = rho_mms(mmsCase, L, x, y);
            V[i][j].u   = uvel_mms(mmsCase, L, x, y);
            V[i][j].v   = vvel_mms(mmsCase, L, x, y);
            V[i][j].P   = press_mms(mmsCase, L, x, y);
        }
    }

    // DO NOT convert here; ghosts are not set yet.
    // (No GlobalPrimitiveToConserved() call in initializeMMS.)
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
  const std::vector<std::vector<double>>& cellVolume,   // <-- add this
  std::vector<std::vector<Conserved>>&   S
) {
  assert(mmsCase==1 || mmsCase==2);
  int Ni = imax + 2*ghost;
  int Nj = jmax + 2*ghost;
  S.assign(Ni, std::vector<Conserved>(Nj, Conserved{0,0,0,0}));

  for (int i = ghost; i < ghost+imax; ++i) {
    for (int j = ghost; j < ghost+jmax; ++j) {
      double x = x_cell[i][j];
      double y = y_cell[i][j];
      double vol = cellVolume[i][j];
      // source density (per volume) * volume -> integral
      S[i][j].rho  = rmassconv (mmsCase, L, x, y) * vol;
      S[i][j].rhou = xmtmconv  (mmsCase, L, x, y) * vol;
      S[i][j].rhov = ymtmconv  (mmsCase, L, x, y) * vol;
      S[i][j].E    = energyconv(mmsCase, gamma, L, x, y) * vol;
    }
  }
}



void applyBoundaryConditions(
    std::vector<std::vector<Conserved>> &U,   // (can be unused if you prefer)
    std::vector<std::vector<Primitive>> &V,
    int mmsCase, double L,
    const std::vector<std::vector<double>> &x_cell,
    const std::vector<std::vector<double>> &y_cell
) {
    int ni = imax + 2 * ghost;
    int nj = jmax + 2 * ghost;

    for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
            // only ghosts
            if (i >= ghost && i < ghost + imax && j >= ghost && j < ghost + jmax) continue;

            double x = x_cell[i][j];
            double y = y_cell[i][j];

            V[i][j].rho = rho_mms(mmsCase, L, x, y);
            V[i][j].u   = uvel_mms(mmsCase, L, x, y);
            V[i][j].v   = vvel_mms(mmsCase, L, x, y);
            V[i][j].P   = press_mms(mmsCase, L, x, y);
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
        for (int i = ghost; i < ghost + imax; ++i) {
          for (int j = ghost; j < ghost + jmax; ++j) {
           const Primitive& Vcell = V[i][j];
            double rho = std::max(Vcell.rho, RHO_MIN);
            double P   = std::max(Vcell.P,   P_MIN);
            double a   = std::sqrt(gamma * P / rho);
            double u = Vcell.u, v = Vcell.v;

            double lambda_iL = std::abs(u * nx_face_i[i][j]   + v * ny_face_i[i][j])   + a;
            double lambda_iR = std::abs(u * nx_face_i[i+1][j] + v * ny_face_i[i+1][j]) + a;
            double lambda_jB = std::abs(u * nx_face_j[i][j]   + v * ny_face_j[i][j])   + a;
            double lambda_jT = std::abs(u * nx_face_j[i][j+1] + v * ny_face_j[i][j+1]) + a;

            double areaSum = lambda_iL * A_face_i[i][j]
                           + lambda_iR * A_face_i[i+1][j]
                           + lambda_jB * A_face_j[i][j]
                           + lambda_jT * A_face_j[i][j+1];

            if (areaSum <= 1e-14) continue; // skip degenerate cell
            double dt_cell = CFL * cellVolume[i][j] / areaSum;
            if (std::isfinite(dt_cell) && dt_cell < dtMin) dtMin = dt_cell;
          }
        }
        return dtMin;
    }


// Area of cell (i,j) from *nodal* corners: (i,j), (i+1,j), (i+1,j+1), (i,j+1)
inline double cellAreaFromNodes(
    int i, int j,
    const std::vector<std::vector<double>>& Xn,
    const std::vector<std::vector<double>>& Yn)
{
    const double x00 = Xn[i  ][j  ], y00 = Yn[i  ][j  ];
    const double x10 = Xn[i+1][j  ], y10 = Yn[i+1][j  ];
    const double x11 = Xn[i+1][j+1], y11 = Yn[i+1][j+1];
    const double x01 = Xn[i  ][j+1], y01 = Yn[i  ][j+1];

    // split quad into two triangles: (00,10,11) + (00,11,01)
    auto tri2 = [](double xa,double ya,double xb,double yb,double xc,double yc){
        return 0.5*std::abs( (xb-xa)*(yc-ya) - (yb-ya)*(xc-xa) );
    };
    return tri2(x00,y00,x10,y10,x11,y11) + tri2(x00,y00,x11,y11,x01,y01);
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
  return (r + std::fabs(r)) / safeDenom(1.0 + r);
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

inline double xi_vanleer(double r){
    return (r + std::fabs(r)) / std::max(1e-12, 1.0 + r);
}
inline double rratio(double up, double um, double small=1e-12){
    double d = up - um; // central jump
    return (std::fabs(d) < small) ? 0.0 : (up - um) / d; // you’ll pass proper args below
}

template<class T> inline void clampPrim(T& W){
    W.rho = std::max(W.rho, RHO_MIN);
    W.P   = std::max(W.P,   P_MIN);
    W.u   = std::clamp(W.u, U_MIN, U_MAX);
    W.v   = std::clamp(W.v, V_MIN, V_MAX);
}

inline void musclI(const std::vector<std::vector<Primitive>>& V, int i,int j,
                   int order,double kappa,bool freeze, Primitive& L, Primitive& R)
{
    const double eps = (order==2 ? 1.0 : 0.0);

    auto rec = [&](auto get, auto set){
        const double LL = get(V[i-1][j]);
        const double C  = get(V[i  ][j]);
        const double R1 = get(V[i+1][j]);
        const double R2 = get(V[i+2][j]);

        const double dC = R1 - C;
        const double dL = C  - LL;
        const double dR = R2 - R1;

        const double rP = (std::fabs(dC) < 1e-12) ? 0.0 : dR / dC;
        const double rM = (std::fabs(dC) < 1e-12) ? 0.0 : dL / dC;

        const double xiP = freeze ? 1.0 : xi_vanleer(rP);
        const double xiM = freeze ? 1.0 : xi_vanleer(rM);

        const double corrL = (eps/4.0) * ((1-kappa)*xiP*(C-LL) + (1+kappa)*xiM*(R1-C));
        const double corrR = (eps/4.0) * ((1-kappa)*xiM*(R2-R1) + (1+kappa)*xiP*(R1-C));

        L.*set = C  + corrL;
        R.*set = R1 - corrR;
    };

    rec([](auto&p){return p.rho;}, &Primitive::rho);
    rec([](auto&p){return p.u;  }, &Primitive::u);
    rec([](auto&p){return p.v;  }, &Primitive::v);
    rec([](auto&p){return p.P;  }, &Primitive::P);

    clampPrim(L);
    clampPrim(R);
}

inline void musclJ(const std::vector<std::vector<Primitive>>& V, int i,int j,
                   int order,double kappa,bool freeze, Primitive& L, Primitive& R)
{
    const double eps = (order==2 ? 1.0 : 0.0);

    auto rec = [&](auto get, auto set){
        const double BB = get(V[i][j-1]);
        const double C  = get(V[i][j  ]);
        const double T1 = get(V[i][j+1]);
        const double T2 = get(V[i][j+2]);

        const double dC = T1 - C;
        const double dB = C  - BB;
        const double dT = T2 - T1;

        const double rP = (std::fabs(dC) < 1e-12) ? 0.0 : dT / dC;
        const double rM = (std::fabs(dC) < 1e-12) ? 0.0 : dB / dC;

        const double xiP = freeze ? 1.0 : xi_vanleer(rP);
        const double xiM = freeze ? 1.0 : xi_vanleer(rM);

        const double corrL = (eps/4.0) * ((1-kappa)*xiP*(C-BB) + (1+kappa)*xiM*(T1-C));
        const double corrR = (eps/4.0) * ((1-kappa)*xiM*(T2-T1) + (1+kappa)*xiP*(T1-C));

        L.*set = C  + corrL;
        R.*set = T1 - corrR;
    };

    rec([](auto&p){return p.rho;}, &Primitive::rho);
    rec([](auto&p){return p.u;  }, &Primitive::u);
    rec([](auto&p){return p.v;  }, &Primitive::v);
    rec([](auto&p){return p.P;  }, &Primitive::P);

    clampPrim(L);
    clampPrim(R);
}

//------------------------------------------------------------------------------
// computeResidualMMS
//
// Builds the residual over each control volume as:
//     R[i][j] = (∑face convective fluxes)  −  S[i][j]
// where S is the *volume-integrated* MMS source, so that the exact MMS solution
// satisfies R = 0.
//
// fluxOrder     : 1 or 2
// kappa         : limiter parameter (e.g. −1, 0, 0.5, 1)
// freezeLimiter : if true, bypasses the limiter (→ pure MUSCL without TVD lim)
// x_cell, y_cell: cell-center coordinates (not used here but kept for parity)
// V             : primitive array (rho,u,v,P) including ghosts
// S             : volume-integrated source terms (same indexing as V)
// R             : output residual, same size as V
//------------------------------------------------------------------------------
void computeResidualMMS(
    int fluxOrder,
    double kappa,
    bool freezeLimiter,
    const std::vector<std::vector<double>>& x_cell,
    const std::vector<std::vector<double>>& y_cell,
    const std::vector<std::vector<Primitive>>& V,
    const std::vector<std::vector<Conserved>>& S,
    std::vector<std::vector<Conserved>>& R
)
{
    int Ni = imax + 2*ghost;
    int Nj = jmax + 2*ghost;

    R.assign(Ni, std::vector<Conserved>(Nj, {0,0,0,0}));

    // Vertical i-faces: include left physical boundary (i = ghost-1)
    for (int i = ghost-1; i <= ghost + imax - 1; ++i) {
        for (int j = ghost; j < ghost + jmax; ++j) {
            Primitive PL, PR;
            musclI(V, i, j, fluxOrder, kappa, freezeLimiter, PL, PR);
            auto F = faceFluxVL2D(PL, PR, nx_face_i[i][j], ny_face_i[i][j], A_face_i[i][j]);
            R[i  ][j] += F;  // left cell
            R[i+1][j] -= F;  // right cell
        }
    }

    // Horizontal j-faces: include bottom physical boundary (j = ghost-1)
    for (int i = ghost; i < ghost + imax; ++i) {
        for (int j = ghost-1; j <= ghost + jmax - 1; ++j) {
            Primitive PL, PR;
            musclJ(V, i, j, fluxOrder, kappa, freezeLimiter, PL, PR);
            auto G = faceFluxVL2D(PL, PR, nx_face_j[i][j], ny_face_j[i][j], A_face_j[i][j]);
            R[i][j  ] += G;  // bottom cell
            R[i][j+1] -= G;  // top cell
        }
    }

    // Subtract the MMS source (volume-integrated) on interior cells
    for (int i = ghost; i < ghost + imax; ++i) {
        for (int j = ghost; j < ghost + jmax; ++j) {
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
    std::vector<std::vector<Conserved>>& R_int,
    double dt,
    const std::vector<std::vector<double>>& cellVolume,
    bool /*debugMode*/, double /*dx*/, double /*dy*/
) {
    const int i0 = ghost, i1 = ghost + imax - 1;
    const int j0 = ghost, j1 = ghost + jmax - 1;

    // Ensure ghosts are consistent at the start of the step
    applyBoundaryConditions(U, V, mmsCase, L, x_cell, y_cell);

    // ===== Stage 1 @ V^n =====
    computeResidualMMS(fluxOrder, kappa, freezeLimiter, x_cell, y_cell, V, S, R_int);
    const auto& R1 = R_int;

    auto U0     = U;
    auto U_star = U0;

    for (int i = i0; i <= i1; ++i)
        for (int j = j0; j <= j1; ++j) {
            double vol = cellVolume[i][j];
            U_star[i][j].rho  = U0[i][j].rho  - (dt/vol) * R1[i][j].rho;
            U_star[i][j].rhou = U0[i][j].rhou - (dt/vol) * R1[i][j].rhou;
            U_star[i][j].rhov = U0[i][j].rhov - (dt/vol) * R1[i][j].rhov;
            U_star[i][j].E    = U0[i][j].E    - (dt/vol) * R1[i][j].E;
        }

    // Apply BCs to U* ghosts (if your BCs write U)
    applyBoundaryConditions(U_star, V /*unused here*/, mmsCase, L, x_cell, y_cell);

    // Build V* from U*
    std::vector<std::vector<Primitive>> V_star(U_star.size(),
                                               std::vector<Primitive>(U_star[0].size()));
    for (int i = 0; i < (int)U_star.size(); ++i)
        for (int j = 0; j < (int)U_star[0].size(); ++j)
            V_star[i][j] = ConservedToPrimitiveCell(U_star[i][j]);

    // 🔧 Critical: apply BCs to V* ghosts as well (Dirichlet primitives for MMS)
    applyBoundaryConditions(U_star /*ignored*/, V_star, mmsCase, L, x_cell, y_cell);

    // ===== Stage 2 @ V* =====
    computeResidualMMS(fluxOrder, kappa, freezeLimiter, x_cell, y_cell, V_star, S, R_int);
    const auto& R2 = R_int;

    // Heun average → U^{n+1}
    for (int i = i0; i <= i1; ++i)
        for (int j = j0; j <= j1; ++j) {
            double vol = cellVolume[i][j];
            U[i][j].rho  = U0[i][j].rho  - 0.5*(dt/vol) * (R1[i][j].rho  + R2[i][j].rho);
            U[i][j].rhou = U0[i][j].rhou - 0.5*(dt/vol) * (R1[i][j].rhou + R2[i][j].rhou);
            U[i][j].rhov = U0[i][j].rhov - 0.5*(dt/vol) * (R1[i][j].rhov + R2[i][j].rhov);
            U[i][j].E    = U0[i][j].E    - 0.5*(dt/vol) * (R1[i][j].E    + R2[i][j].E);
        }

    // Sync primitives from updated U and re-apply BCs for ghosts
    GlobalConservedToPrimitive();                   // V ← U (interior+ghosts)
    applyBoundaryConditions(U, V, mmsCase, L, x_cell, y_cell);
}





struct Residual4 {
    double cont;    // continuity (rho)
    double momx;    // x-momentum (rhou)
    double momy;    // y-momentum (rhov)
    double energy;  // total energy (E)
    double combined; // optional: sqrt(cont^2 + momx^2 + momy^2 + energy^2)
};

// L2 over *interior* cells only
inline Residual4 computeResidualL2_4(const std::vector<std::vector<Conserved>>& R)
{
    double s_rho=0, s_rhou=0, s_rhov=0, s_E=0;
    for (int i = ghost; i < ghost + imax; ++i)
      for (int j = ghost; j < ghost + jmax; ++j) {
        const auto& q = R[i][j];
        s_rho  += q.rho  * q.rho;
        s_rhou += q.rhou * q.rhou;
        s_rhov += q.rhov * q.rhov;
        s_E    += q.E    * q.E;
      }
    Residual4 out;
    out.cont   = std::sqrt(s_rho);
    out.momx   = std::sqrt(s_rhou);
    out.momy   = std::sqrt(s_rhov);
    out.energy = std::sqrt(s_E);
    out.combined = std::sqrt(out.cont*out.cont + out.momx*out.momx
                           + out.momy*out.momy + out.energy*out.energy);
    return out;
}

// Volume-normalize a residual field (R/Vol) cell-wise
inline std::vector<std::vector<Conserved>>
volumeNormalize(const std::vector<std::vector<Conserved>>& Rin,
                const std::vector<std::vector<double>>& cellVolume)
{
    auto Q = Rin;
    for (int i = ghost; i < ghost + imax; ++i)
      for (int j = ghost; j < ghost + jmax; ++j) {
        double vol = cellVolume[i][j];
        Q[i][j].rho  /= vol;
        Q[i][j].rhou /= vol;
        Q[i][j].rhov /= vol;
        Q[i][j].E    /= vol;
      }
    return Q;
}

// Safe component-wise normalization current/base
inline Residual4 normalizeResidual4(const Residual4& cur, const Residual4& base)
{
    auto nd = [](double a, double b){ return a / std::max(b, 1e-300); };
    return {
        nd(cur.cont,   base.cont),
        nd(cur.momx,   base.momx),
        nd(cur.momy,   base.momy),
        nd(cur.energy, base.energy),
        nd(cur.combined, base.combined)
    };
}


//------------------------------------------------------------------------------
// 2D Tecplot dump: writes one POINT‐format zone per call.
//  Assumes globals imax, jmax, ghost, gamma, and arrays
//    x_cell, y_cell  (cell‐center coords, sized Ni×Nj)
//    V                (Primitive), sized Ni×Nj
// are all defined.
//------------------------------------------------------------------------------
// Write nodes for the mesh, and cell-centered primitives (+Mach).
// If writeGhosts=false: writes ONLY physical domain (I=imax+1, J=jmax+1; vars size imax×jmax)
// If writeGhosts=true : writes ghosts too (I=imax+1+2g, J=jmax+1+2g; vars size (imax+2g)×(jmax+2g))
void OutputSolution2D(const std::string &filename, int iter, bool writeGhosts=false) {
    std::ofstream f(filename, std::ios::app);
    if (!f) { std::cerr << "Error: cannot open " << filename << "\n"; return; }

    auto safe = [](double v){ return std::isfinite(v) ? v : -999.9; };
    auto pos  = [](double v, double lo){ return (std::isfinite(v) && v>lo)? v : lo; };

    const int g = ghost;
    int I_nodes = writeGhosts ? (imax+1 + 2*g) : (imax+1);
    int J_nodes = writeGhosts ? (jmax+1 + 2*g) : (jmax+1);
    int I_cells = writeGhosts ? (imax + 2*g)   : (imax);
    int J_cells = writeGhosts ? (jmax + 2*g)   : (jmax);

    // Build nodal array to write (pad nodes only when writing ghosts)
    std::vector<std::vector<double>> Xw, Yw;
    if (writeGhosts) {
        Xw = padWithGhosts2D(x_node, g);
        Yw = padWithGhosts2D(y_node, g);
    }

    if (iter == 0) {
        f << "TITLE = \"2D Euler MMS Solution\"\n";
        f << "VARIABLES = \"X\" \"Y\" \"rho\" \"u\" \"v\" \"P\" \"Mach\"\n";
    }
    f << "ZONE T=\"iter=" << iter << (writeGhosts ? "_ghosts" : "") << "\""
      << ", STRANDID=1, SOLUTIONTIME=" << iter
      << ", I=" << I_nodes << " J=" << J_nodes
      << ", DATAPACKING=BLOCK\n"
      << "VARLOCATION=([3-7]=CELLCENTERED)\n";

    // X nodes
    for (int j=0; j<J_nodes; ++j)
        for (int i=0; i<I_nodes; ++i)
            f << safe(writeGhosts ? Xw[i][j] : x_node[i][j]) << "\n";

    // Y nodes
    for (int j=0; j<J_nodes; ++j)
        for (int i=0; i<I_nodes; ++i)
            f << safe(writeGhosts ? Yw[i][j] : y_node[i][j]) << "\n";

    // Cell-centered: if writing ghosts, dump entire padded V;
    // otherwise write only physical cells (offset by +g,+g)
    for (int j=0; j<J_cells; ++j)
      for (int i=0; i<I_cells; ++i) {
        const auto &W = writeGhosts ? V[i][j] : V[i+g][j+g];
        f << safe(W.rho) << "\n";
      }

    for (int j=0; j<J_cells; ++j)
      for (int i=0; i<I_cells; ++i) {
        const auto &W = writeGhosts ? V[i][j] : V[i+g][j+g];
        f << safe(W.u) << "\n";
      }

    for (int j=0; j<J_cells; ++j)
      for (int i=0; i<I_cells; ++i) {
        const auto &W = writeGhosts ? V[i][j] : V[i+g][j+g];
        f << safe(W.v) << "\n";
      }

    for (int j=0; j<J_cells; ++j)
      for (int i=0; i<I_cells; ++i) {
        const auto &W = writeGhosts ? V[i][j] : V[i+g][j+g];
        f << safe(W.P) << "\n";
      }

    for (int j=0; j<J_cells; ++j)
      for (int i=0; i<I_cells; ++i) {
        const auto &W = writeGhosts ? V[i][j] : V[i+g][j+g];
        double rho = pos(W.rho, RHO_MIN), P = pos(W.P, P_MIN);
        double a = std::sqrt(gamma * P / rho);
        double q = std::sqrt(W.u*W.u + W.v*W.v);
        f << safe(q / (a > 0 ? a : 1e-12)) << "\n";
      }

    f.close();
    std::cout << "[INFO] Wrote zone iter=" << iter
              << (writeGhosts? " (with ghosts)" : " (physical only)") << "\n";
}

inline double specificInternalEnergy(const Conserved& Uc) {
    double rho = std::max(Uc.rho, RHO_MIN);
    double kin = 0.5 * (Uc.rhou*Uc.rhou + Uc.rhov*Uc.rhov) / (rho*rho);
    return Uc.E / rho - kin; // e (per mass)
}

struct FieldStats {
    double minRho, minP, minEint;
    double maxRho, maxP;
    size_t nNaN = 0, nInf = 0;
};

FieldStats check_UV(const char* tag) {
    FieldStats s{+1e99,+1e99,+1e99,-1e99,-1e99,0,0};
    for (int i=ghost; i<ghost+imax; ++i)
      for (int j=ghost; j<ghost+jmax; ++j) {
        const auto& Uc = U[i][j];
        const auto& Vc = V[i][j];
        double e = specificInternalEnergy(Uc);
        auto upd=[&](double v){ if(!std::isfinite(v)) (std::isinf(v)?++s.nInf:++s.nNaN); };
        upd(Uc.rho); upd(Uc.rhou); upd(Uc.rhov); upd(Uc.E);
        upd(Vc.rho); upd(Vc.P);   upd(Vc.u);    upd(Vc.v);

        s.minRho = std::min(s.minRho, Vc.rho);
        s.maxRho = std::max(s.maxRho, Vc.rho);
        s.minP   = std::min(s.minP,   Vc.P);
        s.maxP   = std::max(s.maxP,   Vc.P);
        s.minEint= std::min(s.minEint,e);
      }
    std::cout << "[CHECK] " << tag
              << "  minRho=" << s.minRho
              << "  minP="   << s.minP
              << "  min(e)=" << s.minEint
              << "  NaN/Inf="<< s.nNaN << "/" << s.nInf << "\n";
    return s;
}

static void checkGhostsAgainstMMS(
    const std::vector<std::vector<Primitive>>& V,
    int mmsCase, double L,
    const std::vector<std::vector<double>>& x_cell,
    const std::vector<std::vector<double>>& y_cell,
    const char* where
){
    const int g = ghost;
    auto check = [&](int i, int j, double &eR, double &eU, double &eV, double &eP){
        double x = x_cell[i][j], y = y_cell[i][j];
        Primitive exact{
            rho_mms(mmsCase,L,x,y),
            uvel_mms(mmsCase,L,x,y),
            vvel_mms(mmsCase,L,x,y),
            press_mms(mmsCase,L,x,y)
        };
        eR = std::max(eR, std::abs(V[i][j].rho - exact.rho));
        eU = std::max(eU, std::abs(V[i][j].u   - exact.u));
        eV = std::max(eV, std::abs(V[i][j].v   - exact.v));
        eP = std::max(eP, std::abs(V[i][j].P   - exact.P));
    };

    double eR=0,eU=0,eV=0,eP=0;

    // left/right ghost belts
    for(int q=1; q<=g; ++q){
        int iL = g - q;
        int iR = g + imax - 1 + q;
        for(int j=g; j<g+jmax; ++j){
            check(iL,j,eR,eU,eV,eP);
            check(iR,j,eR,eU,eV,eP);
        }
    }
    // bottom/top ghost belts
    for(int q=1; q<=g; ++q){
        int jB = g - q;
        int jT = g + jmax - 1 + q;
        for(int i=g; i<g+imax; ++i){
            check(i,jB,eR,eU,eV,eP);
            check(i,jT,eR,eU,eV,eP);
        }
    }

    std::cerr << "[BCCHK] " << where
              << "  max|ρ|=" << eR
              << "  |u|="   << eU
              << "  |v|="   << eV
              << "  |P|="   << eP << "\n";
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
    bool debugMode = true;  // set true if you want a simple Cartesian mesh  
    
    double xmin, xmax, ymin, ymax, dx, dy;

    std::vector<std::vector<double>> cellVolume;  // padded volumes (with ghosts)
    
    if (!debugMode) {
     readCurviMeshFromFile(meshFile, x_cell,y_cell,A_face_i,A_face_j,nx_face_i,ny_face_i,nx_face_j,ny_face_j, cellVolume, xmin,xmax,ymin,ymax,dx,dy);


      std::cout << "[INFO] Loaded curvi mesh with imax = " << imax << ", jmax = " << jmax << "\n";
    } else {
    std::ifstream in(meshFile);
    if (!in) { std::cerr << "Error: cannot open " << meshFile << "\n"; std::exit(EXIT_FAILURE); }

    int nz, Ni_nodes, Nj_nodes, kplanes;
    in >> nz >> Ni_nodes >> Nj_nodes >> kplanes;
    // Treat header numbers as NODES; cells = nodes-1
    imax = Ni_nodes - 1;
    jmax = Nj_nodes - 1;

    double Lx = L;  // or set from file/desired size
    double Ly = L;

    buildCartesianDebug(
        imax, jmax, Lx, Ly, ghost,
        x_cell, y_cell,
        A_face_i, A_face_j,
        nx_face_i, ny_face_i,
        nx_face_j, ny_face_j,
        cellVolume,
        x_node, y_node,
        dx, dy
    );

    xmin = 0.0; xmax = Lx; ymin = 0.0; ymax = Ly;

    std::cout << "[INFO] Using Cartesian debug: cells="<<imax<<"x"<<jmax
              << " nodes="<<imax+1<<"x"<<jmax+1
              << " (ghost="<<ghost<<")\n";
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


    // After geometry is built and U/V are resized:
    std::cout << "[INFO] Sample values from x_cell and y_cell:\n";
    int Ishow = std::min(5, (int)x_cell.size());
    int Jshow = std::min(1, (int)(x_cell.empty() ? 0 : x_cell[0].size()));
    for (int i = 0; i < Ishow; ++i) {
        std::cout << "x_cell[" << i << "] = " << x_cell[i][0]
                  << ", y_cell[" << i << "] = " << y_cell[i][0] << "\n";
    }

    // Face lengths
    std::cout << "[INFO] Sample values from A_face_i and A_face_j:\n";
    int Ishow_i = std::min(5, (int)A_face_i.size());
    int Ishow_j = std::min(5, (int)A_face_j.size());
    if (Ishow_i > 0 && !A_face_i[0].empty())
      for (int i = 0; i < Ishow_i; ++i)
        std::cout << "A_face_i[" << i << "] = " << A_face_i[i][0] << "\n";
    if (Ishow_j > 0 && !A_face_j[0].empty())
      for (int i = 0; i < Ishow_j; ++i)
        std::cout << "A_face_j[" << i << "] = " << A_face_j[i][0] << "\n";

    // Normals
    std::cout << "[INFO] Sample values from normal vectors nx/ny:\n";
    int Nxi = std::min(5, (int)nx_face_i.size());
    if (Nxi > 0 && !nx_face_i[0].empty()) {
      for (int i = 0; i < Nxi; ++i) {
       std::cout << "nx_face_i[" << i << "] = " << nx_face_i[i][0]
                  << ", ny_face_i[" << i << "] = " << ny_face_i[i][0] << "\n";
      }
    }
    int Nxj = std::min(5, (int)nx_face_j.size());
    if (Nxj > 0 && !nx_face_j[0].empty()) {
      for (int i = 0; i < Nxj; ++i) {
        std::cout << "nx_face_j[" << i << "] = " << nx_face_j[i][0]
                  << ", ny_face_j[" << i << "] = " << ny_face_j[i][0] << "\n";
      }
    }

    if (!cellVolume.empty() && !cellVolume[0].empty()) {
      std::cout << "[INFO] Sample values from cellVolume array:\n";
      int iBeg = ghost, iEnd = std::min(ghost + imax, (int)cellVolume.size());
      int jBeg = ghost, jEnd = std::min(ghost + jmax, (int)cellVolume[0].size());
      for (int i = iBeg; i < std::min(iBeg+3, iEnd); ++i)
        for (int j = jBeg; j < std::min(jBeg+3, jEnd); ++j)
          std::cout << "cellVolume[" << i << "][" << j << "] = " << cellVolume[i][j] << "\n";
    }

    // Call initializeMMS to fill U/V with the exact primitive → conserved MMS solution:
    int mmsCase = 1;
    initializeMMS(mmsCase, L, x_cell, y_cell);

    // Apply the Dirichlet Boundary Conditions
    applyBoundaryConditions(U, V, mmsCase, L, x_cell, y_cell);

    GlobalPrimitiveToConserved();

    check_UV("after init+BC");

    std::vector<std::vector<Conserved>> S;

    computeSourceTermsMMS(mmsCase, gamma, L, x_cell, y_cell, cellVolume, S);

    GlobalPrimitiveToConserved();

    OutputSolution2D(solFile, 0);

    std::vector<std::vector<Conserved>> R0;
    computeResidualMMS(fluxOrder, kappa, freezeLimiter, x_cell, y_cell, V, S, R0);

    // Baselines (raw and per-volume) for each equation
    Residual4 base_raw4 = computeResidualL2_4(R0);
    Residual4 base_vol4 = computeResidualL2_4(volumeNormalize(R0, cellVolume));

    std::cout << "[INIT] L2 residuals (raw): "
              << "cont="   << base_raw4.cont
              << " momx="  << base_raw4.momx
              << " momy="  << base_raw4.momy
              << " energy="<< base_raw4.energy
              << " | combined=" << base_raw4.combined << "\n";

    std::cout << "[INIT] L2 residuals (perVol): "
              << "cont="   << base_vol4.cont
              << " momx="  << base_vol4.momx
              << " momy="  << base_vol4.momy
              << " energy="<< base_vol4.energy
              << " | combined=" << base_vol4.combined << "\n";


    double dt = computeTimeStep(V, cellVolume, A_face_i, A_face_j, nx_face_i, ny_face_i, nx_face_j, ny_face_j);
    std::cout << "[INFO] Computed time step: dt = " << dt << "\n";
    
    std::vector<std::vector<Conserved>> R_int(imax + 2*ghost,std::vector<Conserved>(jmax + 2*ghost));


        // Pseudo-time marching parameters
    const int maxIter = 100;
    const int writeInterval = 5;
    double tol = 1e-8;             // convergence tolerance
 


    
    for (int iter = 1; iter < maxIter; ++iter) {
    // 1) CFL time step from current V
        double dt = computeTimeStep(
            V, cellVolume,
            A_face_i, A_face_j,
            nx_face_i, ny_face_i, nx_face_j, ny_face_j
        );

    // 2) Advance one RK2 step (updates U and V via your function)
        rungeKutta2Step(
            fluxOrder, kappa, freezeLimiter,
            mmsCase, L,
            S,              // MMS source
            R_int,          // scratch
            dt,
            cellVolume,
            debugMode, dx, dy
        );

        // 3) Residual at updated state
        std::vector<std::vector<Conserved>> R_now;
        computeResidualMMS(fluxOrder, kappa, freezeLimiter, x_cell, y_cell, V, S, R_now);

        // 4) Current L2 norms (raw and per-volume)
        Residual4 cur_raw4 = computeResidualL2_4(R_now);
        Residual4 cur_vol4 = computeResidualL2_4(volumeNormalize(R_now, cellVolume));

        // 5) Normalize component-wise against the *initial* baselines
        Residual4 nraw = normalizeResidual4(cur_raw4, base_raw4);
        Residual4 nvol = normalizeResidual4(cur_vol4, base_vol4);

        // 6) Print ONLY normalized values (use per-volume as your main metric)
        std::cout << "iter=" << iter
                  << "  norm(perVol): "
                  << "cont="   << nvol.cont
                  << " momx="  << nvol.momx
                  << " momy="  << nvol.momy
                  << " energy="<< nvol.energy
                  << " | max=" << std::max(std::max(nvol.cont, nvol.momx),
                                       std::max(nvol.momy, nvol.energy))
                  << "\n";

        // 7) Convergence check on the worst component (safer than combined)
        double conv_metric = std::max(std::max(nvol.cont, nvol.momx),
                                      std::max(nvol.momy, nvol.energy));
        if (conv_metric < tol) {
            std::cout << "Converged at iter=" << iter
                      << " (max normalized per-volume residual=" << conv_metric << ")\n";
            OutputSolution2D(solFile, iter); // final write
            break;
        }

        // 8) Optional: write solution every writeInterval
        if ((iter % writeInterval) == 0) {
            OutputSolution2D(solFile, iter);
        }

    }

    return 0; 
}