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
#include <limits>

using namespace std;
namespace fs = std::filesystem;

#ifndef PI
constexpr double PI = 3.14159265358979323846;
#endif

//-----------------------------------------------------
// Global constants and parameters
//-----------------------------------------------------
constexpr double Rgas  = 287.0;   // [J/(kg·K)]
constexpr double gamma = 1.4;
constexpr double L     = 1.0;     // domain length

static constexpr int ghost = 2;

const double CFL       = 0.01;    // CFL number for time step control
const double epsM      = 0.01;    // Minimum Mach number allowed
const double tolerance = 1e-12;
const double delta     = 1e-6;    // small number used in limiters, etc.

// Nodal coordinates (no ghosts)
std::vector<std::vector<double>> x_node;
std::vector<std::vector<double>> y_node;

static_assert(ghost >= 2, "MUSCL reconstruction needs at least 2 ghost cells.");

//-----------------------------------------------------
// Primitive and Conserved definitions
//-----------------------------------------------------
struct Primitive {
    double rho;  // density
    double u;    // x-velocity
    double v;    // y-velocity
    double P;    // pressure
};

struct Conserved {
    double rho;   // mass
    double rhou;  // x-momentum
    double rhov;  // y-momentum
    double E;     // total energy (per volume)
};

//-----------------------------------------------------
// Global Arrays for the CFD Solver
//-----------------------------------------------------
int imax = 0, jmax = 0;
int kmax; // (unused here)

static std::vector<std::vector<Conserved>> U;
static std::vector<std::vector<Primitive>> V;

// Inline operator overloads for Conserved
inline Conserved operator+(const Conserved& a, const Conserved& b) {
    return { a.rho + b.rho, a.rhou + b.rhou, a.rhov + b.rhov, a.E + b.E };
}
inline Conserved operator-(const Conserved& a, const Conserved& b) {
    return { a.rho - b.rho, a.rhou - b.rhou, a.rhov - b.rhov, a.E - b.E };
}
inline Conserved operator*(double s, const Conserved& b) {
    return { s * b.rho, s * b.rhou, s * b.rhov, s * b.E };
}
inline Conserved operator/(const Conserved& a, double s) {
    return { a.rho / s, a.rhou / s, a.rhov / s, a.E / s };
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
const double U1_MIN = RHO_MIN;
const double U1_MAX = RHO_MAX * 2.0;

const double U2_MIN = RHO_MIN * U_MIN * 2.0;   // rho*u
const double U2_MAX = RHO_MAX * U_MAX * 2.0;

const double U3_MIN = RHO_MIN * V_MIN * 2.0;   // rho*v
const double U3_MAX = RHO_MAX * V_MAX * 2.0;

const double U4_MIN = P_MIN / (gamma - 1.0);  // energy density lower bound
const double U4_MAX = P_MAX / (gamma - 1.0)
                    + 0.5 * RHO_MAX * (U_MAX*U_MAX + V_MAX*V_MAX);

// Clamp helpers
void ApplyLimitsToPrimitive(const std::string& when,
                            std::vector<std::vector<Primitive>>& Varr)
{
    int ni = (int)Varr.size();
    if (ni == 0) return;
    int nj = (int)Varr[0].size();

    std::vector<std::pair<int,int>> clampedCells;
    clampedCells.reserve(5);

    for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
            bool clamped = false;
            auto& cell = Varr[i][j];

            double new_rho = std::clamp(cell.rho, RHO_MIN, RHO_MAX);
            if (new_rho != cell.rho) { clamped = true; cell.rho = new_rho; }

            double new_u = std::clamp(cell.u, U_MIN, U_MAX);
            if (new_u != cell.u) { clamped = true; cell.u = new_u; }

            double new_v = std::clamp(cell.v, V_MIN, V_MAX);
            if (new_v != cell.v) { clamped = true; cell.v = new_v; }

            double new_P = std::clamp(cell.P, P_MIN, P_MAX);
            if (new_P != cell.P) { clamped = true; cell.P = new_P; }

            if (clamped && clampedCells.size() < 5) {
                clampedCells.emplace_back(i, j);
            }
        }
    }

    if (!clampedCells.empty()) {
        std::cerr << "[WARNING] Clamped primitive vars"
                  << (when.empty() ? "" : " (" + when + ")")
                  << " in "
                  << clampedCells.size()
                  << " cell" << (clampedCells.size() > 1 ? "s" : "")
                  << ", e.g. ";
        for (auto [i,j] : clampedCells) std::cerr << "(" << i << "," << j << ") ";
        if (clampedCells.size() == 5) std::cerr << "…";
        std::cerr << "\n";
    }
}

void ApplyLimitsToConserved(const std::string& when,
                            std::vector<std::vector<Conserved>>& Uarr)
{
    int ni = (int)Uarr.size();
    if (ni == 0) return;
    int nj = (int)Uarr[0].size();

    std::vector<std::pair<int,int>> clampedCells;
    clampedCells.reserve(5);

    for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
            bool clamped = false;
            auto& cell = Uarr[i][j];

            double new_rho  = std::clamp(cell.rho,  U1_MIN, U1_MAX);
            if (new_rho != cell.rho) { clamped = true; cell.rho = new_rho; }

            double new_rhou = std::clamp(cell.rhou, U2_MIN, U2_MAX);
            if (new_rhou != cell.rhou) { clamped = true; cell.rhou = new_rhou; }

            double new_rhov = std::clamp(cell.rhov, U3_MIN, U3_MAX);
            if (new_rhov != cell.rhov) { clamped = true; cell.rhov = new_rhov; }

            double new_E    = std::clamp(cell.E,    U4_MIN, U4_MAX);
            if (new_E != cell.E) { clamped = true; cell.E = new_E; }

            if (clamped && clampedCells.size() < 5) {
                clampedCells.emplace_back(i, j);
            }
        }
    }

    if (!clampedCells.empty()) {
        std::cerr << "[WARNING] Clamped conserved vars"
                  << (when.empty() ? "" : " (" + when + ")")
                  << " in "
                  << clampedCells.size()
                  << " cell" << (clampedCells.size() > 1 ? "s" : "")
                  << ", e.g. ";
        for (auto [i,j] : clampedCells) std::cerr << "(" << i << "," << j << ") ";
        if (clampedCells.size() == 5) std::cerr << "…";
        std::cerr << "\n";
    }
}

//-----------------------------------------------------
// Conversion Functions (cell-wise)
//-----------------------------------------------------
Conserved PrimitiveToConserved(const Primitive& Vcell) {
    const double q2   = Vcell.u*Vcell.u + Vcell.v*Vcell.v;
    const double Evol = Vcell.P / (gamma - 1.0) + 0.5 * Vcell.rho * q2;

    Conserved Ucell;
    Ucell.rho  = Vcell.rho;
    Ucell.rhou = Vcell.rho * Vcell.u;
    Ucell.rhov = Vcell.rho * Vcell.v;
    Ucell.E    = Evol;  // total energy per volume
    return Ucell;
}

Primitive ConservedToPrimitiveCell(const Conserved& Ucell) {
    Primitive Vcell;
    Vcell.rho = Ucell.rho;

    if (std::fabs(Ucell.rho) < 1e-12) {
        Vcell.u = 0.0;
        Vcell.v = 0.0;
        // With rho ~ 0, just give a tiny positive pressure
        Vcell.P = P_MIN;
        return Vcell;
    }

    Vcell.u = Ucell.rhou / Ucell.rho;
    Vcell.v = Ucell.rhov / Ucell.rho;
    // U.E = p/(gamma-1) + 0.5*rho*(u^2+v^2)

    const double KE_vol = 0.5 * (Ucell.rhou*Ucell.rhou + Ucell.rhov*Ucell.rhov) / Ucell.rho;
    
    Vcell.P = (gamma - 1.0) * (Ucell.E - 0.5 * KE_vol);

    if (Vcell.P < P_MIN) Vcell.P = P_MIN;
    return Vcell;
}

//-----------------------------------------------------
// Global Conversion Routines
//-----------------------------------------------------
void GlobalConservedToPrimitive() {
    int Ni = imax + 2*ghost, Nj = jmax + 2*ghost;
    V.resize(Ni);
    for (int i = 0; i < Ni; ++i) V[i].resize(Nj);

    for (int i = 0; i < Ni; ++i)
        for (int j = 0; j < Nj; ++j)
            V[i][j] = ConservedToPrimitiveCell(U[i][j]);

    ApplyLimitsToPrimitive("after U→V conversion", V);
}

void GlobalPrimitiveToConserved() {
    int Ni = imax + 2*ghost, Nj = jmax + 2*ghost;
    U.resize(Ni);
    for (int i = 0; i < Ni; ++i) U[i].resize(Nj);

    for (int i = 0; i < Ni; ++i)
        for (int j = 0; j < Nj; ++j)
            U[i][j] = PrimitiveToConserved(V[i][j]);

    ApplyLimitsToConserved("after V→U conversion", U);
}

// ---------- Vector type + operators ----------
struct Vec { double x, y; };

inline Vec operator+(Vec a, Vec b) { return {a.x + b.x, a.y + b.y}; }
inline Vec operator-(Vec a, Vec b) { return {a.x - b.x, a.y - b.y}; }
inline Vec operator*(double s, Vec a) { return {s*a.x, s*a.y}; }
inline Vec operator*(Vec a, double s) { return {s*a.x, s*a.y}; }

inline double dot(Vec a, Vec b){ return a.x*b.x + a.y*b.y; }
inline double norm(Vec a){ return std::hypot(a.x, a.y); }
inline Vec    unit(Vec a){ double n = norm(a); return (n>0) ? (1.0/n)*a : Vec{1,0}; }

// Reflect vector d across a line with tangent t (t can be non‑unit)
inline Vec reflectAcrossTangent(Vec d, Vec t){
    t = unit(t);
    double proj = dot(t, d);
    // d_ref = 2*(t·d)*t - d
    return 2.0*proj*t - d;
}

inline double tri_area(double xa,double ya,double xb,double yb,double xc,double yc){
    return 0.5 * std::abs((xb - xa)*(yc - ya) - (yb - ya)*(xc - xa));
}

static Vec boundaryTangentJ(const std::vector<std::vector<double>>& X,
                            const std::vector<std::vector<double>>& Y,
                            int ib, int j, int jmin, int jmax)
{
    if (j>jmin && j<jmax) return unit({ X[ib][j+1]-X[ib][j-1], Y[ib][j+1]-Y[ib][j-1] });
    if (j==jmin && j+1<=jmax) return unit({ X[ib][j+1]-X[ib][j], Y[ib][j+1]-Y[ib][j] });
    if (j==jmax && j-1>=jmin) return unit({ X[ib][j]-X[ib][j-1], Y[ib][j]-Y[ib][j-1] });
    return {1,0};
}
static Vec boundaryTangentI(const std::vector<std::vector<double>>& X,
                            const std::vector<std::vector<double>>& Y,
                            int i, int jb, int imin, int imax)
{
    if (i>imin && i<imax) return unit({ X[i+1][jb]-X[i-1][jb], Y[i+1][jb]-Y[i-1][jb] });
    if (i==imin && i+1<=imax) return unit({ X[i+1][jb]-X[i][jb], Y[i+1][jb]-Y[i][jb] });
    if (i==imax && i-1>=imin) return unit({ X[i][jb]-X[i-1][jb], Y[i][jb]-Y[i-1][jb] });
    return {1,0};
}

#define GEOMETRIC_MIRROR 0

// Input: interior nodes Xn,Yn of size Ni_nodes x Nj_nodes
// Output: Xg,Yg of size (Ni_nodes+2g) x (Nj_nodes+2g) with ghosts filled
static void buildGhostNodes(
    const std::vector<std::vector<double>>& Xn,
    const std::vector<std::vector<double>>& Yn,
    int imax, int jmax, int g,
    std::vector<std::vector<double>>& Xg,
    std::vector<std::vector<double>>& Yg)
{
    const int Ni = imax+1, Nj = jmax+1;
    Xg.assign(Ni+2*g, std::vector<double>(Nj+2*g));
    Yg.assign(Ni+2*g, std::vector<double>(Nj+2*g));

    // copy interior to center
    for (int j=0;j<Nj;++j)
      for (int i=0;i<Ni;++i){
        Xg[i+g][j+g]=Xn[i][j];
        Yg[i+g][j+g]=Yn[i][j];
      }

#if !GEOMETRIC_MIRROR
    // ---------- Index-space mirror (simple) ----------
    // Left/right belts
    for (int j=g;j<g+Nj;++j)
      for (int q=1;q<=g;++q){
        Xg[g-q][j]       = 2*Xg[g][j]       - Xg[g+q][j];
        Yg[g-q][j]       = 2*Yg[g][j]       - Yg[g+q][j];
        Xg[g+Ni-1+q][j]  = 2*Xg[g+Ni-1][j]  - Xg[g+Ni-1-q][j];
        Yg[g+Ni-1+q][j]  = 2*Yg[g+Ni-1][j]  - Yg[g+Ni-1-q][j];
      }
    // Bottom/top belts
    for (int i=0;i<Ni+2*g;++i)
      for (int q=1;q<=g;++q){
        Xg[i][g-q]       = 2*Xg[i][g]       - Xg[i][g+q];
        Yg[i][g-q]       = 2*Yg[i][g]       - Yg[i][g+q];
        Xg[i][g+Nj-1+q]  = 2*Xg[i][g+Nj-1]  - Xg[i][g+Nj-1-q];
        Yg[i][g+Nj-1+q]  = 2*Yg[i][g+Nj-1]  - Yg[i][g+Nj-1-q];
      }
#else
    // ---------- Geometric mirror (reflect across boundary polyline) ----------
    const int iMin=g, iMax=g+Ni-1, jMin=g, jMax=g+Nj-1;

    // Left boundary (i=iMin)
    for (int j=jMin;j<=jMax;++j){
        Vec Pb{Xg[iMin][j],Yg[iMin][j]};
        Vec t = boundaryTangentJ(Xg,Yg,iMin,j,jMin,jMax);
        for (int k=1;k<=g;++k){
            Vec Pin{Xg[iMin+k][j],Yg[iMin+k][j]};
            Vec d  = Pin - Pb;
            Vec Pg = Pb + reflectAcrossTangent(d, t);
            Xg[iMin-k][j]=Pg.x; Yg[iMin-k][j]=Pg.y;
        }
    }
    // Right boundary (i=iMax)
    for (int j=jMin;j<=jMax;++j){
        Vec Pb{Xg[iMax][j],Yg[iMax][j]};
        Vec t = boundaryTangentJ(Xg,Yg,iMax,j,jMin,jMax);
        for (int k=1;k<=g;++k){
            Vec Pin{Xg[iMax-k][j],Yg[iMax-k][j]};
            Vec d  = Pin - Pb;
            Vec Pg = Pb + reflectAcrossTangent(d, t);
            Xg[iMax+k][j]=Pg.x; Yg[iMax+k][j]=Pg.y;
        }
    }
    // Bottom boundary (j=jMin)
    for (int i=iMin;i<=iMax;++i){
        Vec Pb{Xg[i][jMin],Yg[i][jMin]};
        Vec t = boundaryTangentI(Xg,Yg,i,jMin,iMin,iMax);
        for (int k=1;k<=g;++k){
            Vec Pin{Xg[i][jMin+k],Yg[i][jMin+k]};
            Vec d  = Pin - Pb;
            Vec Pg = Pb + reflectAcrossTangent(d, t);
            Xg[i][jMin-k]=Pg.x; Yg[i][jMin-k]=Pg.y;
        }
    }
    // Top boundary (j=jMax)
    for (int i=iMin;i<=iMax;++i){
        Vec Pb{Xg[i][jMax],Yg[i][jMax]};
        Vec t = boundaryTangentI(Xg,Yg,i,jMax,iMin,iMax);
        for (int k=1;k<=g;++k){
            Vec Pin{Xg[i][jMax-k],Yg[i][jMax-k]};
            Vec d  = Pin - Pb;
            Vec Pg = Pb + reflectAcrossTangent(d, t);
            Xg[i][jMax+k]=Pg.x; Yg[i][jMax+k]=Pg.y;
        }
    }
#endif
}

// Cell centers (imax+2g x jmax+2g)
static void buildCellCentersFromNodesGhost(
    const std::vector<std::vector<double>>& Xg,
    const std::vector<std::vector<double>>& Yg,
    int imax, int jmax, int g,
    std::vector<std::vector<double>>& x_cell,
    std::vector<std::vector<double>>& y_cell)
{
    x_cell.assign(imax+2*g, std::vector<double>(jmax+2*g));
    y_cell.assign(imax+2*g, std::vector<double>(jmax+2*g));
    for (int j=0;j<jmax+2*g;++j)
      for (int i=0;i<imax+2*g;++i){
        double x00=Xg[i][j],     y00=Yg[i][j];
        double x10=Xg[i+1][j],   y10=Yg[i+1][j];
        double x11=Xg[i+1][j+1], y11=Yg[i+1][j+1];
        double x01=Xg[i][j+1],   y01=Yg[i][j+1];
        x_cell[i][j] = 0.25*(x00+x10+x11+x01);
        y_cell[i][j] = 0.25*(y00+y10+y11+y01);
      }
}

// Face lengths ("areas") and outward normals
static void computeFaceGeomFromNodesGhost(
    const std::vector<std::vector<double>>& Xg,
    const std::vector<std::vector<double>>& Yg,
    int imax, int jmax, int g,
    std::vector<std::vector<double>>& A_face_i,
    std::vector<std::vector<double>>& nx_face_i,
    std::vector<std::vector<double>>& ny_face_i,
    std::vector<std::vector<double>>& A_face_j,
    std::vector<std::vector<double>>& nx_face_j,
    std::vector<std::vector<double>>& ny_face_j)
{
    // i-faces: (imax+1+2g) x (jmax+2g)
    A_face_i.assign(imax+1+2*g, std::vector<double>(jmax+2*g));
    nx_face_i.assign(imax+1+2*g, std::vector<double>(jmax+2*g));
    ny_face_i.assign(imax+1+2*g, std::vector<double>(jmax+2*g));
    for (int j=0;j<jmax+2*g;++j)
      for (int i=0;i<imax+1+2*g;++i){
        double dx = Xg[i][j+1]-Xg[i][j];
        double dy = Yg[i][j+1]-Yg[i][j];
        double L  = std::hypot(dx,dy);
        A_face_i[i][j]=L;
        // outward normal w.r.t. "left" cell (i-1): rotate CCW → (dy,-dx)/L
        nx_face_i[i][j] = (L>0)?  dy/L : 0.0;
        ny_face_i[i][j] = (L>0)? -dx/L : 0.0;
      }

    // j-faces: (imax+2g) x (jmax+1+2g)
    A_face_j.assign(imax+2*g, std::vector<double>(jmax+1+2*g));
    nx_face_j.assign(imax+2*g, std::vector<double>(jmax+1+2*g));
    ny_face_j.assign(imax+2*g, std::vector<double>(jmax+1+2*g));
    for (int j=0;j<jmax+1+2*g;++j)
      for (int i=0;i<imax+2*g;++i){
        double dx = Xg[i+1][j]-Xg[i][j];
        double dy = Yg[i+1][j]-Yg[i][j];
        double L  = std::hypot(dx,dy);
        A_face_j[i][j]=L;
        // outward normal w.r.t. "bottom" cell (j-1): rotate CW → (-dy,dx)/L
        nx_face_j[i][j] = (L>0)? -dy/L : 0.0;
        ny_face_j[i][j] = (L>0)?  dx/L : 0.0;
      }
}

// Cell areas (imax+2g x jmax+2g)
static void buildCellAreasFromNodesGhost(
    const std::vector<std::vector<double>>& Xg,
    const std::vector<std::vector<double>>& Yg,
    int imax, int jmax, int g,
    std::vector<std::vector<double>>& cellVol)
{
    cellVol.assign(imax+2*g, std::vector<double>(jmax+2*g));
    for (int j=0;j<jmax+2*g;++j)
      for (int i=0;i<imax+2*g;++i){
        double x00=Xg[i][j],     y00=Yg[i][j];
        double x10=Xg[i+1][j],   y10=Yg[i+1][j];
        double x11=Xg[i+1][j+1], y11=Yg[i+1][j+1];
        double x01=Xg[i][j+1],   y01=Yg[i][j+1];
        cellVol[i][j] = tri_area(x00,y00,x10,y10,x11,y11)
                      + tri_area(x00,y00,x11,y11,x01,y01);
      }
}

void readCurviMeshFromFile(
    const std::string& path,
    std::vector<std::vector<double>>& x_cell_out,   // (imax+2g) x (jmax+2g)
    std::vector<std::vector<double>>& y_cell_out,   // (imax+2g) x (jmax+2g)
    std::vector<std::vector<double>>& A_face_i_out, // (imax+1+2g) x (jmax+2g)
    std::vector<std::vector<double>>& A_face_j_out, // (imax+2g)   x (jmax+1+2g)
    std::vector<std::vector<double>>& nx_face_i_out,
    std::vector<std::vector<double>>& ny_face_i_out,
    std::vector<std::vector<double>>& nx_face_j_out,
    std::vector<std::vector<double>>& ny_face_j_out,
    std::vector<std::vector<double>>& cellVolume_out,
    double& xmin, double& xmax, double& ymin, double& ymax,
    double& dx_approx, double& dy_approx)
{
    std::ifstream in(path);
    if (!in) { std::cerr << "[ERROR] Cannot open mesh: " << path << "\n"; std::exit(EXIT_FAILURE); }

    int nz, Ni_nodes, Nj_nodes, kplanes;
    in >> nz >> Ni_nodes >> Nj_nodes >> kplanes;
    if (!(nz==1 && kplanes>=1)) { std::cerr << "[ERROR] Bad .grd header\n"; std::exit(EXIT_FAILURE); }

    imax = Ni_nodes - 1;
    jmax = Nj_nodes - 1;

    std::vector<std::vector<double>> Xn(Ni_nodes, std::vector<double>(Nj_nodes));
    std::vector<std::vector<double>> Yn(Ni_nodes, std::vector<double>(Nj_nodes));
    double tmp;

    // X planes
    for (int k=0;k<kplanes;++k)
      for (int j=0;j<Nj_nodes;++j)
        for (int i=0;i<Ni_nodes;++i){
            in >> tmp; if (k==0) Xn[i][j]=tmp;
        }
    // Y planes
    for (int k=0;k<kplanes;++k)
      for (int j=0;j<Nj_nodes;++j)
        for (int i=0;i<Ni_nodes;++i){
            in >> tmp; if (k==0) Yn[i][j]=tmp;
        }
    // Z planes (ignored)
    for (int k=0;k<kplanes;++k)
      for (int j=0;j<Nj_nodes;++j)
        for (int i=0;i<Ni_nodes;++i) in >> tmp;

    // bounds (diag)
    xmin=xmax=Xn[0][0]; ymin=ymax=Yn[0][0];
    for (int j=0;j<Nj_nodes;++j)
      for (int i=0;i<Ni_nodes;++i){
        xmin=std::min(xmin,Xn[i][j]); xmax=std::max(xmax,Xn[i][j]);
        ymin=std::min(ymin,Yn[i][j]); ymax=std::max(ymax,Yn[i][j]);
      }
    dx_approx = (imax>0)? (xmax-xmin)/double(imax) : 0.0;
    dy_approx = (jmax>0)? (ymax-ymin)/double(jmax) : 0.0;

    // 1) Build ghosted nodes
    std::vector<std::vector<double>> Xg, Yg;
    buildGhostNodes(Xn,Yn,imax,jmax,ghost,Xg,Yg);

    // (Optional: keep copies if you want to inspect/plot)
    x_node = Xg; y_node = Yg;  // now WITH ghosts

    // 2) Geometry from ghosted nodes
    buildCellCentersFromNodesGhost(Xg,Yg,imax,jmax,ghost, x_cell_out, y_cell_out);
    computeFaceGeomFromNodesGhost(Xg,Yg,imax,jmax,ghost,
                                  A_face_i_out,nx_face_i_out,ny_face_i_out,
                                  A_face_j_out,nx_face_j_out,ny_face_j_out);
    buildCellAreasFromNodesGhost(Xg,Yg,imax,jmax,ghost, cellVolume_out);

    std::cout << "[READ] nodes(int): " << Ni_nodes << " x " << Nj_nodes
              << "  => cells: " << imax << " x " << jmax
              << "  (ghost=" << ghost << ")\n";
}



// New version: ghosted nodes + geometry-from-nodes (consistent with curvi)
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
    std::vector<std::vector<double>>& x_node,  // (imax+1+2g) x (jmax+1+2g)
    std::vector<std::vector<double>>& y_node,  // (imax+1+2g) x (jmax+1+2g)
    double &dx, double &dy
) {
    dx = Lx / static_cast<double>(imax);
    dy = Ly / static_cast<double>(jmax);

    const int Ni_noGhost = imax + 1;
    const int Nj_noGhost = jmax + 1;
    const int Ni = Ni_noGhost + 2*g; // nodes with ghosts
    const int Nj = Nj_noGhost + 2*g;

    // ----- Ghosted nodes -----
    x_node.assign(Ni, std::vector<double>(Nj));
    y_node.assign(Ni, std::vector<double>(Nj));

    // Fill interior node plane at indices [g..g+imax] × [g..g+jmax]
    for (int j=0; j<Nj_noGhost; ++j)
      for (int i=0; i<Ni_noGhost; ++i){
        x_node[i+g][j+g] = i * dx;
        y_node[i+g][j+g] = j * dy;
      }

    // Extrapolate ghost node belts uniformly (index-space mirror is exact here)
    // Left/right
    for (int j=g; j<g+Nj_noGhost; ++j)
      for (int q=1; q<=g; ++q){
        x_node[g-q][j]       = x_node[g][j]       - q*dx;
        y_node[g-q][j]       = y_node[g][j];
        x_node[g+imax+q][j]  = x_node[g+imax][j]  + q*dx;
        y_node[g+imax+q][j]  = y_node[g+imax][j];
      }
    // Bottom/top
    for (int i=0; i<Ni; ++i)
      for (int q=1; q<=g; ++q){
        x_node[i][g-q]       = x_node[i][g];
        y_node[i][g-q]       = y_node[i][g]       - q*dy;
        x_node[i][g+jmax+q]  = x_node[i][g+jmax];
        y_node[i][g+jmax+q]  = y_node[i][g+jmax]  + q*dy;
      }

    // ----- Centers, faces, areas from ghosted nodes -----
    buildCellCentersFromNodesGhost(x_node, y_node, imax, jmax, g, x_cell, y_cell);
    computeFaceGeomFromNodesGhost(   x_node, y_node, imax, jmax, g,
                                     A_face_i, nx_face_i, ny_face_i,
                                     A_face_j, nx_face_j, ny_face_j);
    buildCellAreasFromNodesGhost(    x_node, y_node, imax, jmax, g, cellVolume);

    // (Optional) quick asserts for Cartesian:
    //  - A_face_i ≈ dy, A_face_j ≈ dx
    //  - nx_face_i ≈ +1,0 ; nx_face_j ≈ 0,+1 (note signs are per our convention:
    //    i-face normal points to +i wrt left cell; j-face normal points to +j wrt bottom cell)
}

// ---------------- MMS constants (same layout you used) ----------------
struct MmsParams {
    double rho0, rho_x, rho_y, rho_z, a_rho_x, a_rho_y, a_rho_z;
    double u0,   u_x,   u_y,   u_z,   a_u_x,   a_u_y,   a_u_z;
    double v0,   v_x,   v_y,   v_z,   a_v_x,   a_v_y,   a_v_z;
    double p0,   p_x,   p_y,   p_z,   a_p_x,   a_p_y,   a_p_z;
};

// === MMS parameter sets (definitions) ========================================
const MmsParams mmsSup = {
    /*rho0,rho_x,rho_y,rho_z, a_rho_x,a_rho_y,a_rho_z*/
    1.0,   0.10,  0.08,  0.0,  0,0,0,
    /*u0,  u_x,   u_y,   u_z,  a_u_x, a_u_y, a_u_z*/
    300.0, 20.0,  15.0,  0.0,  0,0,0,
    /*v0,  v_x,   v_y,   v_z,  a_v_x, a_v_y, a_v_z*/
    0.0,   15.0,  12.0,  0.0,  0,0,0,
    /*p0,     p_x,    p_y,    p_z,  a_p_x,a_p_y,a_p_z*/
    101325.0, 8000.0, 6000.0, 0.0,  0,0,0
};

const MmsParams mmsSub = {
    1.0,    0.02,  0.02,  0.0,  0,0,0,
    50.0,   10.0,  8.0,   0.0,  0,0,0,
    0.0,    5.0,   6.0,   0.0,  0,0,0,
    101325.0, 3000.0, 3000.0, 0.0,  0,0,0
};


// ---------------- Analytic MMS fields (match Mathematica) ----------------
inline const MmsParams& Csel(int mmsCase) {
    assert(mmsCase == 1 || mmsCase == 2);
    return (mmsCase == 1 ? mmsSup : mmsSub);
}

inline double rho_mms(int mmsCase, double L, double x, double y) {
    const auto& C = Csel(mmsCase);
    // rho0 + rhox*sin(pi*x/L) + rhoy*cos(pi*y/(2L))
    return C.rho0
         + C.rho_x * std::sin((PI * x) / L)
         + C.rho_y * std::cos((PI * y) / (2.0 * L));
}

inline double uvel_mms(int mmsCase, double L, double x, double y) {
    const auto& C = Csel(mmsCase);
    // u0 + uvelx*sin(3pi*x/(2L)) + uvely*cos(3pi*y/(5L))
    return C.u0
         + C.u_x * std::sin((3.0 * PI * x) / (2.0 * L))
         + C.u_y * std::cos((3.0 * PI * y) / (5.0 * L));
}

inline double vvel_mms(int mmsCase, double L, double x, double y) {
    const auto& C = Csel(mmsCase);
    // v0 + vvelx*cos(pi*x/(2L)) + vvely*sin(2pi*y/(3L))
    return C.v0
         + C.v_x * std::cos((PI * x) / (2.0 * L))
         + C.v_y * std::sin((2.0 * PI * y) / (3.0 * L));
}

inline double press_mms(int mmsCase, double L, double x, double y) {
    const auto& C = Csel(mmsCase);
    // press0 + pressx*cos(2pi*x/L) + pressy*sin(pi*y/L)
    return C.p0
         + C.p_x * std::cos((2.0 * PI * x) / L)
         + C.p_y * std::sin((PI * y) / L);
}

// ---------------- Source terms S(U) (match FortranForm) ----------------
// Continuity (mass)
inline double rmassconv(int mmsCase, double L, double x, double y) {
    const auto& C = Csel(mmsCase);
    const double rhoPiece = C.rho0
        + C.rho_x * std::sin((PI * x) / L)
        + C.rho_y * std::cos((PI * y) / (2.0 * L));
    const double uPiece = C.u0
        + C.u_x * std::sin((3.0 * PI * x) / (2.0 * L))
        + C.u_y * std::cos((3.0 * PI * y) / (5.0 * L));
    const double vPiece = C.v0
        + C.v_x * std::cos((PI * x) / (2.0 * L))
        + C.v_y * std::sin((2.0 * PI * y) / (3.0 * L));

    const double term1 = (3.0 * PI * C.u_x * std::cos((3.0 * PI * x) / (2.0 * L)) * rhoPiece) / (2.0 * L);
    const double term2 = (2.0 * PI * C.v_y * std::cos((2.0 * PI * y) / (3.0 * L)) * rhoPiece) / (3.0 * L);
    const double term3 = (PI * C.rho_x * std::cos((PI * x) / L) * uPiece) / L;
    const double term4 = (PI * C.rho_y * std::sin((PI * y) / (2.0 * L)) * vPiece) / (2.0 * L);
    return term1 + term2 + term3 - term4;
}

// x-momentum
inline double xmtmconv(int mmsCase, double L, double x, double y) {
    const auto& C = Csel(mmsCase);
    const double rhoPiece = C.rho0
        + C.rho_x * std::sin((PI * x) / L)
        + C.rho_y * std::cos((PI * y) / (2.0 * L));
    const double uPiece = C.u0
        + C.u_x * std::sin((3.0 * PI * x) / (2.0 * L))
        + C.u_y * std::cos((3.0 * PI * y) / (5.0 * L));
    const double vPiece = C.v0
        + C.v_x * std::cos((PI * x) / (2.0 * L))
        + C.v_y * std::sin((2.0 * PI * y) / (3.0 * L));

    const double term1 = (3.0 * PI * C.u_x * std::cos((3.0 * PI * x) / (2.0 * L)) * rhoPiece * uPiece) / L;
    const double term2 = (2.0 * PI * C.v_y * std::cos((2.0 * PI * y) / (3.0 * L)) * rhoPiece * uPiece) / (3.0 * L);
    const double term3 = (PI * C.rho_x * std::cos((PI * x) / L) * uPiece * uPiece) / L;
    const double term4 = (2.0 * PI * C.p_x   * std::sin((2.0 * PI * x) / L)) / L; // pressure gradient
    const double term5 = (PI * C.rho_y * uPiece * std::sin((PI * y) / (2.0 * L)) * vPiece) / (2.0 * L);
    const double term6 = (3.0 * PI * C.u_y * rhoPiece * std::sin((3.0 * PI * y) / (5.0 * L)) * vPiece) / (5.0 * L);
    return term1 + term2 + term3 - term4 - term5 - term6;
}

// y-momentum
inline double ymtmconv(int mmsCase, double L, double x, double y) {
    const auto& C = Csel(mmsCase);
    const double rhoPiece = C.rho0
        + C.rho_x * std::sin((PI * x) / L)
        + C.rho_y * std::cos((PI * y) / (2.0 * L));
    const double uPiece = C.u0
        + C.u_x * std::sin((3.0 * PI * x) / (2.0 * L))
        + C.u_y * std::cos((3.0 * PI * y) / (5.0 * L));
    const double vPiece = C.v0
        + C.v_x * std::cos((PI * x) / (2.0 * L))
        + C.v_y * std::sin((2.0 * PI * y) / (3.0 * L));

    const double term1 = (PI * C.p_y * std::cos((PI * y) / L)) / L;
    const double term2 = (PI * C.v_x * std::sin((PI * x) / (2.0 * L)) * rhoPiece * uPiece) / (2.0 * L);
    const double term3 = (3.0 * PI * C.u_x * std::cos((3.0 * PI * x) / (2.0 * L)) * rhoPiece * vPiece) / (2.0 * L);
    const double term4 = (4.0 * PI * C.v_y * std::cos((2.0 * PI * y) / (3.0 * L)) * rhoPiece * vPiece) / (3.0 * L);
    const double term5 = (PI * C.rho_x * std::cos((PI * x) / L) * uPiece * vPiece) / L;
    const double term6 = (PI * C.rho_y * std::sin((PI * y) / (2.0 * L)) * vPiece * vPiece) / (2.0 * L);
    return term1 - term2 + term3 + term4 + term5 - term6;
}

// Energy (per your notebook: E = 0.5(u^2+v^2) + p/[(γ−1)ρ] in the flux)
inline double energyconv(int mmsCase, double gamma, double L, double x, double y) {
    const auto& C = Csel(mmsCase);

    const double rhoPhi = C.rho0
        + C.rho_x * std::sin((PI * x) / L)
        + C.rho_y * std::cos((PI * y) / (2.0 * L));
    const double uPhi = C.u0
        + C.u_x * std::sin((3.0 * PI * x) / (2.0 * L))
        + C.u_y * std::cos((3.0 * PI * y) / (5.0 * L));
    const double vPhi = C.v0
        + C.v_x * std::cos((PI * x) / (2.0 * L))
        + C.v_y * std::sin((2.0 * PI * y) / (3.0 * L));
    const double pPhi = C.p0
        + C.p_x * std::cos((2.0 * PI * x) / L)
        + C.p_y * std::sin((PI * y) / L);

    const double vel2 = uPhi*uPhi + vPhi*vPhi;
    const double Eeq  = 0.5 * vel2 + pPhi / ((gamma - 1.0) * rhoPhi);

    // bracket1: x-coupled terms
    const double A  = -2.0 * PI * C.p_x * std::sin((2.0 * PI * x) / L) / L;
    const double B  = rhoPhi * ( -2.0 * PI * C.p_x * std::sin((2.0 * PI * x) / L) / ((gamma - 1.0) * L * rhoPhi)
                   + ( (3.0 * PI * C.u_x * std::cos((3.0 * PI * x) / (2.0 * L)) * uPhi
                     -       PI * C.v_x * std::sin((PI * x) / (2.0 * L)) * vPhi) / L ) ) * 0.5;
    const double C1 = -PI * C.rho_x * std::cos((PI * x) / L) * pPhi / ((gamma - 1.0) * L * rhoPhi * rhoPhi);
    const double bracket1 = A + B + C1;

    const double term1 = uPhi * bracket1;
    const double term2 = (PI * C.rho_x * std::cos((PI * x) / L) * Eeq) / L;
    const double term3 = (3.0 * PI * C.u_x * std::cos((3.0 * PI * x) / (2.0 * L)) * (pPhi + rhoPhi * Eeq)) / (2.0 * L);
    const double term4 = (2.0 * PI * C.v_y * std::cos((2.0 * PI * y) / (3.0 * L)) * (pPhi + rhoPhi * Eeq)) / (3.0 * L);

    // bracket2: y pressure gradient & rho_y coupling
    const double bracket2 = (PI * C.p_y * std::cos((PI * y) / L)) / L
                          - (PI * C.rho_y * std::sin((PI * y) / (2.0 * L)) * vPhi) / (2.0 * L);
    const double term5 = vPhi * bracket2;

    // bracket3: y-derivatives in energy transport
    const double bracket3 = (PI * C.p_y * std::cos((PI * y) / L)) / ((gamma - 1.0) * L * rhoPhi)
                          + ( (-6.0 * PI * C.u_y * uPhi * std::sin((3.0 * PI * y) / (5.0 * L))) / (5.0 * L)
                              + (4.0 * PI * C.v_y * std::cos((2.0 * PI * y) / (3.0 * L)) * vPhi) / (3.0 * L) ) * 0.5
                          + (PI * C.rho_y * std::sin((PI * y) / (2.0 * L)) * pPhi) / (2.0 * (gamma - 1.0) * L * rhoPhi * rhoPhi);
    const double term6 = rhoPhi * bracket3;

    return term1 + term2 + term3 + term4 + term5 + term6;
}

void init_mms_fs(int mmsCase, double L,
                 const std::vector<std::vector<double>>& x_cell,
                 const std::vector<std::vector<double>>& y_cell)
{
    const auto& C = (mmsCase==1 ? mmsSup : mmsSub);
    const Primitive Wfs{ C.rho0, C.u0, C.v0, C.p0 };

    const int I = imax + 2*ghost, J = jmax + 2*ghost;
    V.assign(I, std::vector<Primitive>(J)); // no freestream yet

    // 1) interior = MMS (cell centers)
    for (int j = ghost; j < J-ghost; ++j)
      for (int i = ghost; i < I-ghost; ++i) {
        const double x = x_cell[i][j], y = y_cell[i][j];
        V[i][j].rho = rho_mms (mmsCase, L, x, y);
        V[i][j].u   = uvel_mms(mmsCase, L, x, y);
        V[i][j].v   = vvel_mms(mmsCase, L, x, y);
        V[i][j].P   = press_mms(mmsCase, L, x, y);
      }

    // 2) ghosts = freestream (all four belts, incl. corners)
    auto setg = [&](int i,int j){ V[i][j] = Wfs; };
    for (int j=0;j<J;++j) for (int q=0;q<ghost;++q){ setg(       q ,j); setg(I-1-q,j); }
    for (int i=0;i<I;++i) for (int q=0;q<ghost;++q){ setg(i,       q ); setg(i, J-1-q); }

    ApplyLimitsToPrimitive("init_mms_fs", V);
    GlobalPrimitiveToConserved();
}

// BC: Dirichlet MMS on ghost belts (evaluate at ghost cell centers)
// call this before computing fluxes / each RK stage
void bc_mms(int mmsCase, double L,
            const std::vector<std::vector<double>>& x_cell,
            const std::vector<std::vector<double>>& y_cell)
{
    const int I = imax + 2*ghost, J = jmax + 2*ghost;

    auto setg = [&](int i, int j){
        const double x = x_cell[i][j], y = y_cell[i][j];
        Primitive W;
        W.rho = rho_mms (mmsCase, L, x, y);
        W.u   = uvel_mms(mmsCase, L, x, y);
        W.v   = vvel_mms(mmsCase, L, x, y);
        W.P   = press_mms(mmsCase, L, x, y);
        V[i][j] = W;
        U[i][j] = PrimitiveToConserved(W); // touch ghosts only
    };

    // left/right belts (incl. corners)
    for (int j=0; j<J; ++j)
      for (int q=0; q<ghost; ++q) { setg(q,j); setg(I-1-q,j); }

    // bottom/top belts
    for (int i=0; i<I; ++i)
      for (int q=0; q<ghost; ++q) { setg(i,q); setg(i,J-1-q); }
}

// dt from CFL over interior cells
double compute_dt_CFL(double CFLnum,
                      const std::vector<std::vector<double>>& A_face_i,
                      const std::vector<std::vector<double>>& nx_face_i,
                      const std::vector<std::vector<double>>& ny_face_i,
                      const std::vector<std::vector<double>>& A_face_j,
                      const std::vector<std::vector<double>>& nx_face_j,
                      const std::vector<std::vector<double>>& ny_face_j,
                      const std::vector<std::vector<double>>& Vol)
{
    const int I = imax + 2*ghost;
    const int J = jmax + 2*ghost;
    const double tiny = 1e-14;

    double dt_min = std::numeric_limits<double>::infinity();

    for (int j = ghost; j < J - ghost; ++j) {
        for (int i = ghost; i < I - ghost; ++i) {
            const Primitive& W = V[i][j];
            const double a = std::sqrt(gamma * W.P / std::max(W.rho, RHO_MIN));

            double sig = 0.0;

            // left i-face (between (i-1,j) and (i,j)) -> face index i
            {
                const double nx = nx_face_i[i][j], ny = ny_face_i[i][j];
                const double Af = A_face_i[i][j];
                const double un = W.u * nx + W.v * ny;
                sig += (std::fabs(un) + a) * Af;
            }
            // right i-face -> face index i+1
            {
                const double nx = nx_face_i[i+1][j], ny = ny_face_i[i+1][j];
                const double Af = A_face_i[i+1][j];
                const double un = W.u * nx + W.v * ny;
                sig += (std::fabs(un) + a) * Af;
            }
            // bottom j-face (between (i,j-1) and (i,j)) -> face index j
            {
                const double nx = nx_face_j[i][j], ny = ny_face_j[i][j];
                const double Af = A_face_j[i][j];
                const double un = W.u * nx + W.v * ny;
                sig += (std::fabs(un) + a) * Af;
            }
            // top j-face -> face index j+1
            {
                const double nx = nx_face_j[i][j+1], ny = ny_face_j[i][j+1];
                const double Af = A_face_j[i][j+1];
                const double un = W.u * nx + W.v * ny;
                sig += (std::fabs(un) + a) * Af;
            }

            const double vol = std::max(Vol[i][j], tiny);
            const double dt_cell = CFLnum * vol / std::max(sig, tiny);
            if (dt_cell < dt_min) dt_min = dt_cell;
        }
    }

    // fallback guard
    if (!std::isfinite(dt_min) || dt_min <= 0.0) {
        dt_min = 1e-8; // very small safe fallback
    }
    return dt_min;
}

// =========================== 2D van Leer flux (with MUSCL) ===========================

// --- split coefficients (van Leer 1982) ---
inline double Cplus (double M){ if (M<=-1) return 0.0; else if (M<1) return 0.25*(M+1)*(M+1); else return M; }
inline double Cminus(double M){ if (M<=-1) return M;   else if (M<1) return -0.25*(M-1)*(M-1); else return 0.0; }
inline double Dplus (double M){ if (M<=-1) return 0.0; else if (M<1) return 0.25*(M+1)*(M+1)*(-M+2.0); else return 1.0; }
inline double Dminus(double M){ if (M<=-1) return 1.0; else if (M<1) return -0.25*(M-1)*(M-1)*(-M-2.0); else return 0.0; }

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

// ---- small helper: clamp reconstructed primitives for safety ----
inline void clamp_primitive(Primitive& W){
    if (W.rho < RHO_MIN) W.rho = RHO_MIN;
    if (W.P   < P_MIN  ) W.P   = P_MIN;
    // (optional) clamp velocities to your U/V bounds if you like
}

// van Leer normal flux from L/R primitive states at a face with normal (nx,ny)
inline Conserved vl_flux_n(const Primitive& WL, const Primitive& WR, double nx, double ny)
{
    // left state
    const double aL  = std::sqrt(gamma * WL.P / WL.rho);
    const double unL = WL.u*nx + WL.v*ny;
    const double ML  = unL / aL;
    const double CpL = Cplus(ML), DpL = Dplus(ML);
    const double qL2 = WL.u*WL.u + WL.v*WL.v;
    const double HL  = (gamma/(gamma-1.0))*(WL.P/WL.rho) + 0.5*qL2;

    // right state
    const double aR  = std::sqrt(gamma * WR.P / WR.rho);
    const double unR = WR.u*nx + WR.v*ny;
    const double MR  = unR / aR;
    const double CmR = Cminus(MR), DmR = Dminus(MR);
    const double qR2 = WR.u*WR.u + WR.v*WR.v;
    const double HR  = (gamma/(gamma-1.0))*(WR.P/WR.rho) + 0.5*qR2;  // <- fixed

    // convective mass fluxes (per unit length)
    const double mL = WL.rho * aL * CpL;
    const double mR = WR.rho * aR * CmR;

    Conserved Fn;
    Fn.rho  = mL + mR;
    Fn.rhou = mL*WL.u + mR*WR.u + (DpL*WL.P + DmR*WR.P)*nx;
    Fn.rhov = mL*WL.v + mR*WR.v + (DpL*WL.P + DmR*WR.P)*ny;
    Fn.E    = mL*HL   + mR*HR;
    return Fn;
}


// ---- MUSCL (κ-scheme) reconstruction on an i-face (between cells i-1 and i) ----
inline void muscl_reconstruct_i(const std::vector<std::vector<Primitive>>& V,
                                int i_face, int j,
                                int order, double kappa, bool freezeLimiter,
                                Primitive& WL, Primitive& WR)
{
    const int I = (int)V.size();
    const int il = i_face - 1;   // left cell
    const int ir = i_face;       // right cell

    auto recon_scalar = [&](auto getter)->std::pair<double,double>{
        const double v_i   = getter(il  ,j);
        const double v_ip1 = getter(ir  ,j);
        double vL = v_i, vR = v_ip1;

        if (order == 2 && il-2>=0 && ir+1<I){
            const double v_im1 = getter(il-1,j);
            const double v_ip2 = getter(ir+1,j);

            const double d0   = safeDenom(v_ip1 - v_i);
            const double r_p  = (v_ip2 - v_ip1) / d0;    // r+ at i+1/2
            const double r_m  = (v_i   - v_im1) / d0;    // r- at i+1/2
            const double xi_p = freezeLimiter ? 1.0 : xi(r_p);
            const double xi_m = freezeLimiter ? 1.0 : xi(r_m);

            const double d1   = safeDenom(v_ip2 - v_ip1);
            const double r_m_ip1  = (v_ip1 - v_i) / d1;
            const double xi_m_ip1 = freezeLimiter ? 1.0 : xi(r_m_ip1);

            const double eps = 1.0;
            vL = v_i   + (eps/4.0) * ( (1.0-kappa)*xi_p*(v_i   - v_im1) + (1.0+kappa)*xi_m*(v_ip1 - v_i  ) );
            vR = v_ip1 - (eps/4.0) * ( (1.0-kappa)*xi_m_ip1*(v_ip2 - v_ip1) + (1.0+kappa)*xi_p*(v_ip1 - v_i) );
        }
        return {vL, vR};
    };


    auto get_rho = [&](int ii,int jj){ return V[ii][jj].rho; };
    auto get_u   = [&](int ii,int jj){ return V[ii][jj].u;   };
    auto get_v   = [&](int ii,int jj){ return V[ii][jj].v;   };
    auto get_P   = [&](int ii,int jj){ return V[ii][jj].P;   };

    auto [rhoL,rhoR] = recon_scalar(get_rho);
    auto [uL,uR]     = recon_scalar(get_u);
    auto [vL,vR]     = recon_scalar(get_v);
    auto [pL,pR]     = recon_scalar(get_P);

    WL = {rhoL, uL, vL, pL}; clamp_primitive(WL);
    WR = {rhoR, uR, vR, pR}; clamp_primitive(WR);
}


inline void muscl_reconstruct_j(const std::vector<std::vector<Primitive>>& V,
                                int i, int j_face,
                                int order, double kappa, bool freezeLimiter,
                                Primitive& WB, Primitive& WT)
{
    const int J = (int)V[0].size();
    const int jb = j_face - 1;   // bottom cell
    const int jt = j_face;       // top cell

    auto recon_scalar = [&](auto getter)->std::pair<double,double>{
        const double v_j   = getter(i, jb);
        const double v_jp1 = getter(i, jt);

        if (order == 2 && jb-1 >= 0 && jt+1 < J) {
            const double v_jm1 = getter(i, jb-1);
            const double v_jp2 = getter(i, jt+1);

            const double d0   = safeDenom(v_jp1 - v_j);
            const double r_p  = (v_jp2 - v_jp1) / d0;      // r+ at j+1/2
            const double r_m  = (v_j   - v_jm1) / d0;      // r- at j+1/2
            const double xi_p = freezeLimiter ? 1.0 : xi(r_p);
            const double xi_m = freezeLimiter ? 1.0 : xi(r_m);

            const double d1       = safeDenom(v_jp2 - v_jp1);
            const double r_m_jp1  = (v_jp1 - v_j) / d1;    // r- at j+3/2
            const double xi_m_jp1 = freezeLimiter ? 1.0 : xi(r_m_jp1);

            const double eps = 1.0;
            double vB = v_j   + (eps/4.0)*((1.0-kappa)*xi_p*(v_j   - v_jm1) +
                                           (1.0+kappa)*xi_m*(v_jp1 - v_j));
            double vT = v_jp1 - (eps/4.0)*((1.0-kappa)*xi_m_jp1*(v_jp2 - v_jp1) +
                                           (1.0+kappa)*xi_p     *(v_jp1 - v_j));
            return {vB, vT};
        }
        return std::pair<double,double>{v_j, v_jp1};
    };

    auto get_rho = [&](int ii,int jj){ return V[ii][jj].rho; };
    auto get_u   = [&](int ii,int jj){ return V[ii][jj].u;   };
    auto get_v   = [&](int ii,int jj){ return V[ii][jj].v;   };
    auto get_P   = [&](int ii,int jj){ return V[ii][jj].P;   };

    auto [rhoB,rhoT] = recon_scalar(get_rho);
    auto [uB,uT]     = recon_scalar(get_u);
    auto [vBv,vTv]   = recon_scalar(get_v);
    auto [pB,pT]     = recon_scalar(get_P);

    WB = {rhoB, uB, vBv, pB}; clamp_primitive(WB);
    WT = {rhoT, uT, vTv, pT}; clamp_primitive(WT);
}


// ---- Compute normal fluxes (per unit length) on all faces ----
// Fi: size (imax+1+2g) x (jmax+2g)    (i-faces)
// Fj: size (imax+2g)   x (jmax+1+2g)  (j-faces)
// We only fill usable interior/ghost-interior faces: i=1..I-1, j=0..J-1 and j=1..J-1, i=0..I-1.
void compute_fluxes_vl(const std::vector<std::vector<Primitive>>& V,
                       const std::vector<std::vector<double>>& nx_face_i,
                       const std::vector<std::vector<double>>& ny_face_i,
                       const std::vector<std::vector<double>>& nx_face_j,
                       const std::vector<std::vector<double>>& ny_face_j,
                       int order, double kappa, bool freezeLimiter,
                       std::vector<std::vector<Conserved>>& Fi,
                       std::vector<std::vector<Conserved>>& Fj)
{
    const int I = imax + 2*ghost;
    const int J = jmax + 2*ghost;

    Fi.assign(I+1, std::vector<Conserved>(J, {0,0,0,0}));
    Fj.assign(I,   std::vector<Conserved>(J+1, {0,0,0,0}));

    // i-faces: between cells (i-1,j) and (i,j); use MUSCL along i
    for (int j=0; j<J; ++j){
        for (int i=1; i<=I-1; ++i){
            Primitive WL, WR;
            // fall back to 1st-order at very outer bands to avoid out-of-range
            const bool ok2 = (order==2 && i-2>=0 && i+1<I);
            if (ok2) muscl_reconstruct_i(V, i, j, 2, kappa, freezeLimiter, WL, WR);
            else     { WL = V[i-1][j]; WR = V[i][j]; clamp_primitive(WL); clamp_primitive(WR); }

            const double nx = nx_face_i[i][j], ny = ny_face_i[i][j];
            Fi[i][j] = vl_flux_n(WL, WR, nx, ny); // per unit length
        }
    }

    // j-faces: between cells (i,j-1) and (i,j); use MUSCL along j
    for (int j=1; j<=J-1; ++j){
        for (int i=0; i<I; ++i){
            Primitive WB, WT;
            const bool ok2 = (order==2 && j-2>=0 && j+1<J);
            if (ok2) muscl_reconstruct_j(V, i, j, 2, kappa, freezeLimiter, WB, WT);
            else     { WB = V[i][j-1]; WT = V[i][j]; clamp_primitive(WB); clamp_primitive(WT); }

            const double nx = nx_face_j[i][j], ny = ny_face_j[i][j];
            Fj[i][j] = vl_flux_n(WB, WT, nx, ny); // per unit length
        }
    }
}

// 1) flux divergence -> R (Fi,Fj are per-unit-length)
void add_flux_divergence_to_R(
    const std::vector<std::vector<Conserved>>& Fi, // (I+1) x J
    const std::vector<std::vector<Conserved>>& Fj, // I x (J+1)
    const std::vector<std::vector<double>>& A_face_i,
    const std::vector<std::vector<double>>& A_face_j,
    std::vector<std::vector<Conserved>>& R)
{
    const int I = imax + 2*ghost, J = jmax + 2*ghost;

    // i-faces: normal is outward from left cell (i-1) toward right cell (i)
    for (int j=0; j<J; ++j){
        for (int i=1; i<=I-1; ++i){
            const double Af = A_face_i[i][j];
            const Conserved Fn = Fi[i][j];

            // left cell gets +Fn*Af, right cell gets −Fn*Af
            R[i-1][j].rho  += Fn.rho  * Af;
            R[i-1][j].rhou += Fn.rhou * Af;
            R[i-1][j].rhov += Fn.rhov * Af;
            R[i-1][j].E    += Fn.E    * Af;

            R[i  ][j].rho  -= Fn.rho  * Af;
            R[i  ][j].rhou -= Fn.rhou * Af;
            R[i  ][j].rhov -= Fn.rhov * Af;
            R[i  ][j].E    -= Fn.E    * Af;
        }
    }

    // j-faces: normal is outward from bottom cell (j-1) toward top cell (j)
    for (int j=1; j<=J-1; ++j){
        for (int i=0; i<I; ++i){
            const double Af = A_face_j[i][j];
            const Conserved Fn = Fj[i][j];

            // bottom cell gets +Fn*Af, top cell gets −Fn*Af
            R[i][j-1].rho  += Fn.rho  * Af;
            R[i][j-1].rhou += Fn.rhou * Af;
            R[i][j-1].rhov += Fn.rhov * Af;
            R[i][j-1].E    += Fn.E    * Af;

            R[i][j  ].rho  -= Fn.rho  * Af;
            R[i][j  ].rhou -= Fn.rhou * Af;
            R[i][j  ].rhov -= Fn.rhov * Af;
            R[i][j  ].E    -= Fn.E    * Af;
        }
    }
}


// 2) add MMS source at cell centers (interior only)
// 2) add MMS source at cell centers (interior only)
void add_mms_source_to_R(int mmsCase, double L,
                         const std::vector<std::vector<double>>& x_cell,
                         const std::vector<std::vector<double>>& y_cell,
                         const std::vector<std::vector<double>>& Vol,
                         std::vector<std::vector<Conserved>>& R)
{
    const int I = imax + 2*ghost, J = jmax + 2*ghost;
    for (int j=ghost; j<J-ghost; ++j){
        for (int i=ghost; i<I-ghost; ++i){
            const double x = x_cell[i][j], y = y_cell[i][j];
            // NOTE: subtract S * Vol (was +=)
            R[i][j].rho  += rmassconv (mmsCase, L, x, y) * Vol[i][j];
            R[i][j].rhou += xmtmconv  (mmsCase, L, x, y) * Vol[i][j];
            R[i][j].rhov += ymtmconv  (mmsCase, L, x, y) * Vol[i][j];
            R[i][j].E    += energyconv(mmsCase, gamma, L, x, y) * Vol[i][j];
        }
    }
}


// Build G(U) = R/Vol   where  U_t = -G(U).
// Uses: bc_mms, compute_fluxes_vl, add_flux_divergence_to_R, add_mms_source_to_R.
void build_rhs(int mmsCase, double L,
               const std::vector<std::vector<double>>& A_face_i,
               const std::vector<std::vector<double>>& nx_face_i,
               const std::vector<std::vector<double>>& ny_face_i,
               const std::vector<std::vector<double>>& A_face_j,
               const std::vector<std::vector<double>>& nx_face_j,
               const std::vector<std::vector<double>>& ny_face_j,
               const std::vector<std::vector<double>>& x_cell,
               const std::vector<std::vector<double>>& y_cell,
               const std::vector<std::vector<double>>& Vol,
               int order, double kappa, bool freezeLimiter,
               std::vector<std::vector<Conserved>>& G)   // out: same size as U/V
{
    const int I = imax + 2*ghost, J = jmax + 2*ghost;

    // 1) Dirichlet MMS on ghost cells (states needed for face fluxes)
    bc_mms(mmsCase, L, x_cell, y_cell);

    // 2) per-unit-length face fluxes
    std::vector<std::vector<Conserved>> Fi, Fj;
    compute_fluxes_vl(V, nx_face_i, ny_face_i, nx_face_j, ny_face_j,
                      order, kappa, freezeLimiter, Fi, Fj);

    // 3) residual from fluxes × areas
    std::vector<std::vector<Conserved>> R(I, std::vector<Conserved>(J, {0,0,0,0}));
    add_flux_divergence_to_R(Fi, Fj, A_face_i, A_face_j, R);

    // 4) add MMS sources at cell centers (interior only)
    add_mms_source_to_R(mmsCase, L, x_cell, y_cell, Vol, R);

    // 5) G = R/Vol  (we don’t touch ghosts)
    G.assign(I, std::vector<Conserved>(J, {0,0,0,0}));
    for (int j=ghost; j<J-ghost; ++j){
        for (int i=ghost; i<I-ghost; ++i){
            G[i][j].rho  = R[i][j].rho  / Vol[i][j];
            G[i][j].rhou = R[i][j].rhou / Vol[i][j];
            G[i][j].rhov = R[i][j].rhov / Vol[i][j];
            G[i][j].E    = R[i][j].E    / Vol[i][j];
        }
    }
}

// Returns the dt it used (you pass dt in).
double rk2_step(double dt, int mmsCase, double L,
                const std::vector<std::vector<double>>& A_face_i,
                const std::vector<std::vector<double>>& nx_face_i,
                const std::vector<std::vector<double>>& ny_face_i,
                const std::vector<std::vector<double>>& A_face_j,
                const std::vector<std::vector<double>>& nx_face_j,
                const std::vector<std::vector<double>>& ny_face_j,
                const std::vector<std::vector<double>>& x_cell,
                const std::vector<std::vector<double>>& y_cell,
                const std::vector<std::vector<double>>& Vol,
                int order, double kappa, bool freezeLimiter)
{
    const int I = imax + 2*ghost, J = jmax + 2*ghost;

    // Save U^n
    auto U0 = U;

    // k1 = G(U^n)
    std::vector<std::vector<Conserved>> k1;
    build_rhs(mmsCase, L, A_face_i, nx_face_i, ny_face_i,
              A_face_j, nx_face_j, ny_face_j, x_cell, y_cell, Vol,
              order, kappa, freezeLimiter, k1);

    // Stage 1: U1 = U0 - dt * k1
    for (int j=ghost; j<J-ghost; ++j){
        for (int i=ghost; i<I-ghost; ++i){
            U[i][j].rho  = U0[i][j].rho  - dt * k1[i][j].rho;
            U[i][j].rhou = U0[i][j].rhou - dt * k1[i][j].rhou;
            U[i][j].rhov = U0[i][j].rhov - dt * k1[i][j].rhov;
            U[i][j].E    = U0[i][j].E    - dt * k1[i][j].E;
        }
    }
    GlobalConservedToPrimitive();

    // k2 = G(U1)
    std::vector<std::vector<Conserved>> k2;
    build_rhs(mmsCase, L, A_face_i, nx_face_i, ny_face_i,
              A_face_j, nx_face_j, ny_face_j, x_cell, y_cell, Vol,
              order, kappa, freezeLimiter, k2);

    // Final: U^{n+1} = 0.5*U0 + 0.5*(U1 - dt*k2)
    for (int j=ghost; j<J-ghost; ++j){
        for (int i=ghost; i<I-ghost; ++i){
            Conserved tmp;
            tmp.rho  = U[i][j].rho  - dt * k2[i][j].rho;
            tmp.rhou = U[i][j].rhou - dt * k2[i][j].rhou;
            tmp.rhov = U[i][j].rhov - dt * k2[i][j].rhov;
            tmp.E    = U[i][j].E    - dt * k2[i][j].E;

            U[i][j].rho  = 0.5*(U0[i][j].rho  + tmp.rho );
            U[i][j].rhou = 0.5*(U0[i][j].rhou + tmp.rhou);
            U[i][j].rhov = 0.5*(U0[i][j].rhov + tmp.rhov);
            U[i][j].E    = 0.5*(U0[i][j].E    + tmp.E   );
        }
    }
    GlobalConservedToPrimitive();
    ApplyLimitsToConserved("rk2_step", U);  // optional safety
    ApplyLimitsToPrimitive("rk2_step", V);

    return dt;
}

double rk4_step(double dt, int mmsCase, double L,
                const std::vector<std::vector<double>>& A_face_i,
                const std::vector<std::vector<double>>& nx_face_i,
                const std::vector<std::vector<double>>& ny_face_i,
                const std::vector<std::vector<double>>& A_face_j,
                const std::vector<std::vector<double>>& nx_face_j,
                const std::vector<std::vector<double>>& ny_face_j,
                const std::vector<std::vector<double>>& x_cell,
                const std::vector<std::vector<double>>& y_cell,
                const std::vector<std::vector<double>>& Vol,
                int order, double kappa, bool freezeLimiter)
{
    const int I = imax + 2*ghost, J = jmax + 2*ghost;

    // Save U^n
    auto U0 = U;

    // k1
    std::vector<std::vector<Conserved>> k1, k2, k3, k4;

    build_rhs(mmsCase, L, A_face_i, nx_face_i, ny_face_i,
              A_face_j, nx_face_j, ny_face_j, x_cell, y_cell, Vol,
              order, kappa, freezeLimiter, k1);

    // Utmp = U0 - 0.5*dt*k1
    U = U0;
    for (int j=ghost; j<J-ghost; ++j){
        for (int i=ghost; i<I-ghost; ++i){
            U[i][j].rho  -= 0.5*dt * k1[i][j].rho;
            U[i][j].rhou -= 0.5*dt * k1[i][j].rhou;
            U[i][j].rhov -= 0.5*dt * k1[i][j].rhov;
            U[i][j].E    -= 0.5*dt * k1[i][j].E;
        }
    }
    GlobalConservedToPrimitive();

    // k2
    build_rhs(mmsCase, L, A_face_i, nx_face_i, ny_face_i,
              A_face_j, nx_face_j, ny_face_j, x_cell, y_cell, Vol,
              order, kappa, freezeLimiter, k2);

    // Utmp = U0 - 0.5*dt*k2
    U = U0;
    for (int j=ghost; j<J-ghost; ++j){
        for (int i=ghost; i<I-ghost; ++i){
            U[i][j].rho  -= 0.5*dt * k2[i][j].rho;
            U[i][j].rhou -= 0.5*dt * k2[i][j].rhou;
            U[i][j].rhov -= 0.5*dt * k2[i][j].rhov;
            U[i][j].E    -= 0.5*dt * k2[i][j].E;
        }
    }
    GlobalConservedToPrimitive();

    // k3
    build_rhs(mmsCase, L, A_face_i, nx_face_i, ny_face_i,
              A_face_j, nx_face_j, ny_face_j, x_cell, y_cell, Vol,
              order, kappa, freezeLimiter, k3);

    // Utmp = U0 - dt*k3
    U = U0;
    for (int j=ghost; j<J-ghost; ++j){
        for (int i=ghost; i<I-ghost; ++i){
            U[i][j].rho  -= dt * k3[i][j].rho;
            U[i][j].rhou -= dt * k3[i][j].rhou;
            U[i][j].rhov -= dt * k3[i][j].rhov;
            U[i][j].E    -= dt * k3[i][j].E;
        }
    }
    GlobalConservedToPrimitive();

    // k4
    build_rhs(mmsCase, L, A_face_i, nx_face_i, ny_face_i,
              A_face_j, nx_face_j, ny_face_j, x_cell, y_cell, Vol,
              order, kappa, freezeLimiter, k4);

    // Combine: U^{n+1} = U0 - dt*(k1 + 2k2 + 2k3 + k4)/6
    U = U0;
    for (int j=ghost; j<J-ghost; ++j){
        for (int i=ghost; i<I-ghost; ++i){
            const double s1 = k1[i][j].rho  + 2.0*k2[i][j].rho  + 2.0*k3[i][j].rho  + k4[i][j].rho;
            const double s2 = k1[i][j].rhou + 2.0*k2[i][j].rhou + 2.0*k3[i][j].rhou + k4[i][j].rhou;
            const double s3 = k1[i][j].rhov + 2.0*k2[i][j].rhov + 2.0*k3[i][j].rhov + k4[i][j].rhov;
            const double s4 = k1[i][j].E    + 2.0*k2[i][j].E    + 2.0*k3[i][j].E    + k4[i][j].E;

            U[i][j].rho  -= (dt/6.0) * s1;
            U[i][j].rhou -= (dt/6.0) * s2;
            U[i][j].rhov -= (dt/6.0) * s3;
            U[i][j].E    -= (dt/6.0) * s4;
        }
    }
    GlobalConservedToPrimitive();
    ApplyLimitsToConserved("rk4_step", U);  // optional safety clamps
    ApplyLimitsToPrimitive("rk4_step", V);

    return dt;
}

// compute RMS (L2) and Linf per component over INTERIOR cells only
struct CompNorms {
    double l2_cont=0, l2_momx=0, l2_momy=0, l2_energy=0;
    double linf_cont=0, linf_momx=0, linf_momy=0, linf_energy=0;
};

inline CompNorms compute_component_norms(const std::vector<std::vector<Conserved>>& A)
{
    const int I = imax + 2*ghost, J = jmax + 2*ghost;
    double s1=0, s2=0, s3=0, s4=0; // sum of squares
    double m1=0, m2=0, m3=0, m4=0; // max abs
    std::size_t N = 0;

    for (int j=ghost; j<J-ghost; ++j){
        for (int i=ghost; i<I-ghost; ++i){
            const auto& a = A[i][j];
            s1 += a.rho  * a.rho;   m1 = std::max(m1, std::fabs(a.rho ));
            s2 += a.rhou * a.rhou;  m2 = std::max(m2, std::fabs(a.rhou));
            s3 += a.rhov * a.rhov;  m3 = std::max(m3, std::fabs(a.rhov));
            s4 += a.E    * a.E;     m4 = std::max(m4, std::fabs(a.E   ));
            ++N;
        }
    }
    const double invN = (N>0)? 1.0/double(N) : 0.0;

    CompNorms n;
    n.l2_cont   = std::sqrt(s1*invN);
    n.l2_momx   = std::sqrt(s2*invN);
    n.l2_momy   = std::sqrt(s3*invN);
    n.l2_energy = std::sqrt(s4*invN);
    n.linf_cont   = m1;
    n.linf_momx   = m2;
    n.linf_momy   = m3;
    n.linf_energy = m4;
    return n;
}

// pretty printer (single line). Set reset_base=true on first call to get relative norms.
inline void print_residuals(const std::vector<std::vector<Conserved>>& G,
                            int iter, bool reset_base=false,
                            const char* tag="RHS")
{
    static bool have_base=false;
    static CompNorms base;
    const CompNorms cur = compute_component_norms(G);

    if (reset_base || !have_base){ base=cur; have_base=true; }

    auto rel = [](double v, double b){ const double tiny=1e-30; return v/std::max(b,tiny); };

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "[IT " << std::setw(6) << iter << "] " << tag
              << "  L2: cont=" << cur.l2_cont
              << "  momx=" << cur.l2_momx
              << "  momy=" << cur.l2_momy
              << "  energy=" << cur.l2_energy
              << "   |  Linf: cont=" << cur.linf_cont
              << "  momx=" << cur.linf_momx
              << "  momy=" << cur.linf_momy
              << "  energy=" << cur.linf_energy
              << "   |  Rel(L2): "
              << "cont=" << rel(cur.l2_cont, base.l2_cont)
              << "  momx=" << rel(cur.l2_momx, base.l2_momx)
              << "  momy=" << rel(cur.l2_momy, base.l2_momy)
              << "  energy=" << rel(cur.l2_energy, base.l2_energy)
              << "\n";
}

// Writes a Tecplot BLOCK zone. If writeGhosts=true, it writes the padded grid
// (nodes with ghosts) + all cell-centered fields including ghosts.
// Adds a cell-centered "ghost" mask: 1 in ghost belt, 0 in physical interior.
void OutputSolution2D(const std::string &filename, int iter, bool writeGhosts=false)
{
    std::ofstream f(filename, std::ios::app);
    if (!f) { std::cerr << "Error: cannot open " << filename << "\n"; return; }

    auto safe = [](double v){ return std::isfinite(v) ? v : -999.9; };
    auto pos  = [](double v, double lo){ return (std::isfinite(v) && v>lo)? v : lo; };

    const int g = ghost;

    // --- Build nodal arrays to write ---
    // If we’re writing ghosts, construct ghost nodes via your geometric routine.
    // Otherwise, write the interior nodal grid directly.
    std::vector<std::vector<double>> Xw, Yw;
    if (writeGhosts) {
        // uses your function from earlier messages
        buildGhostNodes(x_node, y_node, imax, jmax, g, Xw, Yw);
    } else {
        Xw = x_node; Yw = y_node;
    }

    const int I_nodes = (int)Xw.size();
    const int J_nodes = I_nodes ? (int)Xw[0].size() : 0;
    const int I_cells = I_nodes - 1;          // Tecplot structured rule
    const int J_cells = J_nodes - 1;

    // --- Header (BLOCK, X/Y nodal, rest cell-centered) ---
    if (iter == 0) {
        f << "TITLE = \"2D Euler MMS Solution\"\n";
        // Note: include a ghost mask var at the end
        f << "VARIABLES = \"X\" \"Y\" \"rho\" \"u\" \"v\" \"P\" \"Mach\" \"ghost\"\n";
    }
    f << "ZONE T=\"iter=" << iter << (writeGhosts ? "_ghosts" : "_phys")
      << "\", STRANDID=1, SOLUTIONTIME=" << iter
      << ", I=" << I_nodes << " J=" << J_nodes
      << ", DATAPACKING=BLOCK\n"
      << "VARLOCATION=([3-8]=CELLCENTERED)\n";

    // --- Blocks: X (nodal), then Y (nodal) ---
    for (int j=0; j<J_nodes; ++j)
        for (int i=0; i<I_nodes; ++i)
            f << safe(Xw[i][j]) << "\n";

    for (int j=0; j<J_nodes; ++j)
        for (int i=0; i<I_nodes; ++i)
            f << safe(Yw[i][j]) << "\n";

    // --- Cell-centered blocks: rho, u, v, P, Mach, ghost ---
    // When NOT writing ghosts, take values from interior slice V[i+g][j+g].
    // When writing ghosts, dump the whole padded V as-is.
    for (int j=0; j<J_cells; ++j)
        for (int i=0; i<I_cells; ++i) {
            const Primitive &W = writeGhosts ? V[i][j] : V[i+g][j+g];
            f << safe(W.rho) << "\n";
        }

    for (int j=0; j<J_cells; ++j)
        for (int i=0; i<I_cells; ++i) {
            const Primitive &W = writeGhosts ? V[i][j] : V[i+g][j+g];
            f << safe(W.u) << "\n";
        }

    for (int j=0; j<J_cells; ++j)
        for (int i=0; i<I_cells; ++i) {
            const Primitive &W = writeGhosts ? V[i][j] : V[i+g][j+g];
            f << safe(W.v) << "\n";
        }

    for (int j=0; j<J_cells; ++j)
        for (int i=0; i<I_cells; ++i) {
            const Primitive &W = writeGhosts ? V[i][j] : V[i+g][j+g];
            f << safe(W.P) << "\n";
        }

    for (int j=0; j<J_cells; ++j)
        for (int i=0; i<I_cells; ++i) {
            const Primitive &W = writeGhosts ? V[i][j] : V[i+g][j+g];
            double rho = pos(W.rho, RHO_MIN), P = pos(W.P, P_MIN);
            double a   = std::sqrt(gamma * P / rho);
            double q   = std::sqrt(W.u*W.u + W.v*W.v);
            f << safe(q / (a > 0 ? a : 1e-12)) << "\n";
        }

    // ghost mask (cell-centered): 1 = ghost, 0 = interior
    for (int j=0; j<J_cells; ++j)
        for (int i=0; i<I_cells; ++i) {
            int gi = i, gj = j; // indices into padded arrays when writeGhosts=true
            if (!writeGhosts) { gi = i + g; gj = j + g; }
            const int I0 = g, I1 = g + imax - 1;
            const int J0 = g, J1 = g + jmax - 1;
            const int isGhost = (gi < I0 || gi > I1 || gj < J0 || gj > J1) ? 1 : 0;
            f << isGhost << "\n";
        }

    f.close();
    std::cout << "[INFO] Wrote Tecplot zone iter=" << iter
              << (writeGhosts? " (with ghosts)" : " (physical only)") << "\n";
}

int main() {
    // --------------- run controls ---------------
    const int    fluxOrder     = 1;      // 1 or 2 (MUSCL)
    const double kappa         = 0.0;    // MUSCL κ in [-1,1]
    const bool   freezeLimiter = false;  // true => no limiting (ξ=1)
    const bool   useRK4        = false;  // false->RK2, true->RK4
    const int    maxIter       = 1000;
    const int    writeEvery    = 50;     // Tecplot frequency

    // --------------- output file ---------------
    std::string outFolder = "OutputFiles";
    if (!fs::exists(outFolder)) {
        fs::create_directory(outFolder);
        std::cout << "[INFO] Created folder: " << outFolder << "\n";
    }
    const std::string solFile = outFolder + "/MMS_Solution.dat";
    { std::ofstream(solFile, std::ios::out).close(); } // truncate

    // --------------- mesh / geometry ---------------
    const std::string meshFile =
        R"(G:\Shared drives\AOE Lab7\Monica Shanmugam\MS\CFD Proj Grids\Grids\curviliniear-grids\curv2d9.grd)";
    const bool debugMode = false;

    // geometry containers
    std::vector<std::vector<double>>
        x_cell, y_cell,                     // (imax+2g) x (jmax+2g) cell centers
        A_face_i, A_face_j,                 // lengths on i- and j-faces
        nx_face_i, ny_face_i,               // unit normals on i-faces
        nx_face_j, ny_face_j,               // unit normals on j-faces
        cellVolume;                         // cell areas (with ghosts)

    double xmin=0, xmax=0, ymin=0, ymax=0, dx=0, dy=0;

    if (!debugMode) {
        readCurviMeshFromFile(meshFile,
            x_cell, y_cell,
            A_face_i, A_face_j,
            nx_face_i, ny_face_i,
            nx_face_j, ny_face_j,
            cellVolume,
            xmin, xmax, ymin, ymax,
            dx, dy);
        std::cout << "[INFO] Curvilinear mesh loaded: imax="<<imax<<" jmax="<<jmax
                  << " (ghost="<<ghost<<")\n";
    } else {
        // for debug: use dims from the .grd header so sizes stay consistent
        std::ifstream in(meshFile);
        if (!in) { std::cerr << "Error: cannot open " << meshFile << "\n"; return EXIT_FAILURE; }
        int nz, Ni_nodes, Nj_nodes, kplanes;
        in >> nz >> Ni_nodes >> Nj_nodes >> kplanes;
        imax = Ni_nodes - 1;
        jmax = Nj_nodes - 1;

        const double Lx = L, Ly = L;
        buildCartesianDebug(imax, jmax, Lx, Ly, ghost,
            x_cell, y_cell,
            A_face_i, A_face_j,
            nx_face_i, ny_face_i,
            nx_face_j, ny_face_j,
            cellVolume,
            x_node, y_node,   // will be filled with GHOSTED nodes
            dx, dy);
        xmin=0; xmax=Lx; ymin=0; ymax=Ly;
        std::cout << "[INFO] Cartesian debug: cells="<<imax<<"x"<<jmax
                  << " nodes="<<imax+1<<"x"<<jmax+1
                  << " (ghost="<<ghost<<")\n";
    }

    // sizes
    const int I = imax + 2*ghost;
    const int J = jmax + 2*ghost;

    // --------------- init (MMS interior; freestream ghosts) ---------------
    const int mmsCase = 1;
    init_mms_fs(mmsCase, L, x_cell, y_cell); // fills V (interior MMS, ghosts FS) + U

    // --------------- prepare nodes for Tecplot writer ---------------
    // OutputSolution2D expects interior nodes when writeGhosts=false and
    // will rebuild ghosts itself when writeGhosts=true. De-ghost x_node/y_node:
    {
        std::vector<std::vector<double>> Xint(imax+1, std::vector<double>(jmax+1));
        std::vector<std::vector<double>> Yint(imax+1, std::vector<double>(jmax+1));
        for (int j=0; j<=jmax; ++j)
            for (int i=0; i<=imax; ++i) {
                Xint[i][j] = x_node[i+ghost][j+ghost];
                Yint[i][j] = y_node[i+ghost][j+ghost];
            }
        x_node.swap(Xint);
        y_node.swap(Yint);
    }

    // --------------- baseline RHS + write initial zones ---------------
    // baseline RHS (with MMS Dirichlet ghosts)
    std::vector<std::vector<Conserved>> G0;
    build_rhs(mmsCase, L,
              A_face_i, nx_face_i, ny_face_i,
              A_face_j, nx_face_j, ny_face_j,
              x_cell, y_cell, cellVolume,
              /*order=*/fluxOrder, /*kappa=*/kappa, /*freezeLimiter=*/freezeLimiter,
              G0);
    print_residuals(G0, /*iter=*/0, /*reset_base=*/true, "RHS");

    // write both zones at iter=0
    OutputSolution2D(solFile, /*iter=*/0, /*writeGhosts=*/false);
    OutputSolution2D(solFile, /*iter=*/0, /*writeGhosts=*/true);

    // --------------- time loop ---------------
    for (int iter = 1; iter <= maxIter; ++iter) {
        // CFL dt
        const double dt = compute_dt_CFL(CFL,
                                         A_face_i, nx_face_i, ny_face_i,
                                         A_face_j, nx_face_j, ny_face_j,
                                         cellVolume);

        // advance
        if (useRK4) {
            rk4_step(dt, mmsCase, L,
                     A_face_i, nx_face_i, ny_face_i,
                     A_face_j, nx_face_j, ny_face_j,
                     x_cell, y_cell, cellVolume,
                     /*order=*/fluxOrder, /*kappa=*/kappa, /*freezeLimiter=*/freezeLimiter);
        } else {
            rk2_step(dt, mmsCase, L,
                     A_face_i, nx_face_i, ny_face_i,
                     A_face_j, nx_face_j, ny_face_j,
                     x_cell, y_cell, cellVolume,
                     /*order=*/fluxOrder, /*kappa=*/kappa, /*freezeLimiter=*/freezeLimiter);
        }

        // current RHS and residual print
        std::vector<std::vector<Conserved>> Gnow;
        build_rhs(mmsCase, L,
                  A_face_i, nx_face_i, ny_face_i,
                  A_face_j, nx_face_j, ny_face_j,
                  x_cell, y_cell, cellVolume,
                  /*order=*/fluxOrder, /*kappa=*/kappa, /*freezeLimiter=*/freezeLimiter,
                  Gnow);
        print_residuals(Gnow, iter, /*reset_base=*/false, "RHS");

        // write solution periodically
        if (iter % writeEvery == 0) {
            OutputSolution2D(solFile, iter, /*writeGhosts=*/false);
            OutputSolution2D(solFile, iter, /*writeGhosts=*/true);
        }
    }

    std::cout << "[DONE]\n";
    return 0;
}