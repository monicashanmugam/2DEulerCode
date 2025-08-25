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

// ---------------- Geometry arrays (padded with ghosts) ----------------
std::vector<std::vector<double>> x_cell;
std::vector<std::vector<double>> y_cell;
std::vector<std::vector<double>> A_face_i; // (imax+1+2g) x (jmax+2g)
std::vector<std::vector<double>> A_face_j; // (imax+2g)   x (jmax+1+2g)
std::vector<std::vector<double>> nx_face_i, ny_face_i;
std::vector<std::vector<double>> nx_face_j, ny_face_j;

struct ErrorData2D {
    double dx;   // (approx.) grid spacing in i-direction
    double dy;   // (approx.) grid spacing in j-direction
    double eP;   // L2 error for pressure
    double eRho; // L2 error for density
    double eU;   // L2 error for u
    double eV;   // L2 error for v
};

// ---------------- Small helper ----------------
inline double tri_area(double xa,double ya,double xb,double yb,double xc,double yc){
    return 0.5 * std::abs((xb - xa)*(yc - ya) - (yb - ya)*(xc - xa));
}

// ---------------- Optional: compute faces/normals from nodal X,Y ----------------
// Shapes (no ghosts): X,Y are (imax+1) x (jmax+1)
void computeMeshGeometry(
    const std::vector<std::vector<double>>& X,
    const std::vector<std::vector<double>>& Y,
    std::vector<std::vector<double>>& A_face_i_out, // (imax+1) x jmax
    std::vector<std::vector<double>>& A_face_j_out, // imax x (jmax+1)
    std::vector<std::vector<double>>& nx_face_i_out,
    std::vector<std::vector<double>>& ny_face_i_out,
    std::vector<std::vector<double>>& nx_face_j_out,
    std::vector<std::vector<double>>& ny_face_j_out
){
    A_face_i_out.assign(imax+1, std::vector<double>(jmax));
    nx_face_i_out.assign(imax+1, std::vector<double>(jmax));
    ny_face_i_out.assign(imax+1, std::vector<double>(jmax));

    for (int j = 0; j < jmax; ++j) {
        for (int i = 0; i <= imax; ++i) {
            const double dx = X[i][j+1] - X[i][j];
            const double dy = Y[i][j+1] - Y[i][j];
            const double L  = std::hypot(dx, dy);

            A_face_i_out[i][j] = L;
            // outward normal to the *left* cell: rotate (dx,dy) CCW 90° → ( dy, -dx )
            nx_face_i_out[i][j] =  dy / L;
            ny_face_i_out[i][j] = -dx / L;
        }
    }

    A_face_j_out.assign(imax, std::vector<double>(jmax+1));
    nx_face_j_out.assign(imax, std::vector<double>(jmax+1));
    ny_face_j_out.assign(imax, std::vector<double>(jmax+1));

    for (int j = 0; j <= jmax; ++j) {
        for (int i = 0; i < imax; ++i) {
            const double dx = X[i+1][j] - X[i][j];
            const double dy = Y[i+1][j] - Y[i][j];
            const double L  = std::hypot(dx, dy);

            A_face_j_out[i][j] = L;
            // outward normal to the *bottom* cell: rotate (dx,dy) CW 90° → ( -dy, dx )
            nx_face_j_out[i][j] = -dy / L;
            ny_face_j_out[i][j] =  dx / L;
        }
    }
}

// ---------------- Curvilinear reader + padding (ghosts) ----------------
// Input .grd header:
//   nz
//   Ni_nodes  Nj_nodes  kplanes
// Then X plane(s), Y plane(s), Z plane(s) (Z ignored). We use plane 0.
void readCurviMeshFromFile(
    const std::string& path,
    std::vector<std::vector<double>>& x_cell_out,
    std::vector<std::vector<double>>& y_cell_out,
    std::vector<std::vector<double>>& A_face_i_out,
    std::vector<std::vector<double>>& A_face_j_out,
    std::vector<std::vector<double>>& nx_face_i_out,
    std::vector<std::vector<double>>& ny_face_i_out,
    std::vector<std::vector<double>>& nx_face_j_out,
    std::vector<std::vector<double>>& ny_face_j_out,
    std::vector<std::vector<double>>& cellVolume_out, // (imax+2g) x (jmax+2g) after padding
    double& xmin, double& xmax, double& ymin, double& ymax,
    double& dx_approx, double& dy_approx
){
    std::ifstream in(path);
    if (!in) {
        std::cerr << "[ERROR] Cannot open mesh file: " << path << "\n";
        std::exit(EXIT_FAILURE);
    }

    int nz, Ni_nodes, Nj_nodes, kplanes;
    in >> nz >> Ni_nodes >> Nj_nodes >> kplanes;
    if (!(nz == 1 && kplanes >= 1)) {
        std::cerr << "[ERROR] Bad .grd header (nz or kplanes)\n";
        std::exit(EXIT_FAILURE);
    }

    imax = Ni_nodes - 1;
    jmax = Nj_nodes - 1;

    // Read nodal X, Y (take plane 0)
    std::vector<std::vector<double>> Xn(Ni_nodes, std::vector<double>(Nj_nodes));
    std::vector<std::vector<double>> Yn(Ni_nodes, std::vector<double>(Nj_nodes));

    double tmp;
    for (int k = 0; k < kplanes; ++k)
        for (int j = 0; j < Nj_nodes; ++j)
            for (int i = 0; i < Ni_nodes; ++i) {
                in >> tmp;
                if (k == 0) Xn[i][j] = tmp;
            }

    for (int k = 0; k < kplanes; ++k)
        for (int j = 0; j < Nj_nodes; ++j)
            for (int i = 0; i < Ni_nodes; ++i) {
                in >> tmp;
                if (k == 0) Yn[i][j] = tmp;
            }

    // Skip Z planes
    for (int k = 0; k < kplanes; ++k)
        for (int j = 0; j < Nj_nodes; ++j)
            for (int i = 0; i < Ni_nodes; ++i) in >> tmp;

    // Save nodal arrays (no ghosts) for output/diagnostics
    x_node = Xn;
    y_node = Yn;

    // Bounding box + crude dx,dy estimates (for diagnostics only)
    xmin = xmax = Xn[0][0];
    ymin = ymax = Yn[0][0];
    for (int j = 0; j < Nj_nodes; ++j) {
        for (int i = 0; i < Ni_nodes; ++i) {
            xmin = std::min(xmin, Xn[i][j]);
            xmax = std::max(xmax, Xn[i][j]);
            ymin = std::min(ymin, Yn[i][j]);
            ymax = std::max(ymax, Yn[i][j]);
        }
    }

    // ----- cell centers (interior, no ghosts) -----
    std::vector<std::vector<double>> xc(imax, std::vector<double>(jmax));
    std::vector<std::vector<double>> yc(imax, std::vector<double>(jmax));
    for (int j = 0; j < jmax; ++j)
        for (int i = 0; i < imax; ++i) {
            xc[i][j] = 0.25*(Xn[i][j] + Xn[i+1][j] + Xn[i][j+1] + Xn[i+1][j+1]);
            yc[i][j] = 0.25*(Yn[i][j] + Yn[i+1][j] + Yn[i][j+1] + Yn[i+1][j+1]);
        }

    // ----- interior face lengths & outward normals (no ghosts) -----
    std::vector<std::vector<double>> Afi(imax+1, std::vector<double>(jmax));
    std::vector<std::vector<double>> nxi(imax+1, std::vector<double>(jmax));
    std::vector<std::vector<double>> nyi(imax+1, std::vector<double>(jmax));
    for (int j = 0; j < jmax; ++j) {
        for (int i = 0; i <= imax; ++i) {
            const double dx = Xn[i][j+1] - Xn[i][j];
            const double dy = Yn[i][j+1] - Yn[i][j];
            const double facelength  = std::hypot(dx, dy);
            Afi[i][j] = facelength;
            nxi[i][j] =  dy / facelength; // outward from left cell
            nyi[i][j] = -dx / facelength;
        }
    }

    std::vector<std::vector<double>> Afj(imax, std::vector<double>(jmax+1));
    std::vector<std::vector<double>> nxj(imax, std::vector<double>(jmax+1));
    std::vector<std::vector<double>> nyj(imax, std::vector<double>(jmax+1));
    for (int j = 0; j <= jmax; ++j) {
        for (int i = 0; i < imax; ++i) {
            const double dx = Xn[i+1][j] - Xn[i][j];
            const double dy = Yn[i+1][j] - Yn[i][j];
            const double facelength  = std::hypot(dx, dy);
            Afj[i][j] = facelength;
            nxj[i][j] = -dy / facelength; // outward from bottom cell
            nyj[i][j] =  dx / facelength;
        }
    }

    // ----- interior cell areas from nodal corners (no ghosts) -----
    std::vector<std::vector<double>> vol(imax, std::vector<double>(jmax));
    for (int j = 0; j < jmax; ++j) {
        for (int i = 0; i < imax; ++i) {
            const double x00 = Xn[i  ][j  ], y00 = Yn[i  ][j  ];
            const double x10 = Xn[i+1][j  ], y10 = Yn[i+1][j  ];
            const double x11 = Xn[i+1][j+1], y11 = Yn[i+1][j+1];
            const double x01 = Xn[i  ][j+1], y01 = Yn[i  ][j+1];
            vol[i][j] = tri_area(x00,y00,x10,y10,x11,y11)
                      + tri_area(x00,y00,x11,y11,x01,y01);
        }
    }

    // ----- padding helpers (mirror), then flip normals in ghosts -----
    auto padScalar = [&](const std::vector<std::vector<double>>& in)->std::vector<std::vector<double>> {
        const int g = ghost;
        if (in.empty() || in[0].empty() || g <= 0) return in;
        const int I = (int)in.size();
        const int J = (int)in[0].size();
        std::vector<std::vector<double>> out(I+2*g, std::vector<double>(J+2*g));
        // copy interior
        for (int i = 0; i < I; ++i)
            for (int j = 0; j < J; ++j)
                out[i+g][j+g] = in[i][j];
        // left/right belts (mirror)
        for (int j = g; j < g+J; ++j)
            for (int q = 0; q < g; ++q) {
                out[g-1-q][j]    = out[g+q][j];
                out[g+I+q][j]    = out[g+I-1-q][j];
            }
        // bottom/top belts (mirror)
        for (int i = 0; i < I+2*g; ++i)
            for (int q = 0; q < g; ++q) {
                out[i][g-1-q]    = out[i][g+q];
                out[i][g+J+q]    = out[i][g+J-1-q];
            }
        return out;
    };

    auto padNormals = [&](std::vector<std::vector<double>>& nx,
                          std::vector<std::vector<double>>& ny) {
        const int g = ghost;
        if (nx.empty() || nx[0].empty() || g <= 0) return;
        const int I = (int)nx.size();
        const int J = (int)nx[0].size();
        // left/right belts: mirror and flip sign
        for (int j = g; j <= J-g-1; ++j) {
            for (int q = 0; q < g; ++q) {
                int il = 2*g - 1 - q;
                int ir = I - g - 1 - q;
                nx[q][j]        = -nx[il][j]; ny[q][j]        = -ny[il][j];
                nx[I-g+q][j]    = -nx[ir][j]; ny[I-g+q][j]    = -ny[ir][j];
            }
        }
        // bottom/top belts: mirror and flip sign
        for (int i = 0; i < I; ++i) {
            for (int q = 0; q < g; ++q) {
                int jb = 2*g - 1 - q;
                int jt = J - g - 1 - q;
                nx[i][q]        = -nx[i][jb]; ny[i][q]        = -ny[i][jb];
                nx[i][J-g+q]    = -nx[i][jt]; ny[i][J-g+q]    = -ny[i][jt];
            }
        }
    };

    // ----- pad all arrays used by the solver -----
    x_cell_out     = padScalar(xc);
    y_cell_out     = padScalar(yc);
    A_face_i_out   = padScalar(Afi);
    nx_face_i_out  = padScalar(nxi);
    ny_face_i_out  = padScalar(nyi);
    A_face_j_out   = padScalar(Afj);
    nx_face_j_out  = padScalar(nxj);
    ny_face_j_out  = padScalar(nyj);
    cellVolume_out = padScalar(vol);

    // flip normals in ghost belts so “outward-from-cell” convention stays correct
    padNormals(nx_face_i_out, ny_face_i_out);
    padNormals(nx_face_j_out, ny_face_j_out);

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
    // ----- Nodes (no ghosts) -----
    x_node.assign(imax+1, std::vector<double>(jmax+1));
    y_node.assign(imax+1, std::vector<double>(jmax+1));

    dx = Lx / static_cast<double>(imax);
    dy = Ly / static_cast<double>(jmax);

    for (int j = 0; j <= jmax; ++j) {
        for (int i = 0; i <= imax; ++i) {
            x_node[i][j] = i * dx;
            y_node[i][j] = j * dy;
        }
    }

    // ----- Padded sizes -----
    const int Ni = imax + 2*g;
    const int Nj = jmax + 2*g;

    // ----- Cell centers (padded) -----
    x_cell.assign(Ni, std::vector<double>(Nj));
    y_cell.assign(Ni, std::vector<double>(Nj));

    for (int j = 0; j < Nj; ++j) {
        for (int i = 0; i < Ni; ++i) {
            // physical interior is i∈[g..g+imax-1], j∈[g..g+jmax-1]
            const double xc = (static_cast<double>(i) - g + 0.5) * dx;
            const double yc = (static_cast<double>(j) - g + 0.5) * dy;
            x_cell[i][j] = xc;
            y_cell[i][j] = yc;
        }
    }

    // ----- Faces (padded) — constant lengths, axis-aligned normals -----
    // i-faces: (Ni+1) × Nj
    A_face_i.assign(Ni+1, std::vector<double>(Nj, dy));
    nx_face_i.assign(Ni+1, std::vector<double>(Nj, 1.0));
    ny_face_i.assign(Ni+1, std::vector<double>(Nj, 0.0));

    // j-faces: Ni × (Nj+1)
    A_face_j.assign(Ni, std::vector<double>(Nj+1, dx));
    nx_face_j.assign(Ni, std::vector<double>(Nj+1, 0.0));
    ny_face_j.assign(Ni, std::vector<double>(Nj+1, 1.0));

    // ----- Cell volumes (padded) -----
    cellVolume.assign(Ni, std::vector<double>(Nj, dx * dy));
}






