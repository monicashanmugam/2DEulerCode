#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <cassert>
#include <iomanip>

using namespace std;
namespace fs = std::filesystem;

//-----------------------------------------------------
// Global constants and parameters
//-----------------------------------------------------
constexpr double L      = 1.0;  // Domain length
static constexpr int ghost = 2;
int imax, jmax;

//-----------------------------------------------------
// Global mesh arrays
//-----------------------------------------------------
vector<vector<double>> x_cell, y_cell;
vector<vector<double>> A_face_i, A_face_j;
vector<vector<double>> nx_face_i, ny_face_i;
vector<vector<double>> nx_face_j, ny_face_j;

//-----------------------------------------------------
// Variable definitions
//-----------------------------------------------------
struct Primitive { double rho, u, v, P; };
struct Conserved { double rho, rhou, rhov, E; };

vector<vector<Primitive>> V;
vector<vector<Conserved>> U;

//-----------------------------------------------------
// Conversion: Primitive -> Conserved
//-----------------------------------------------------
inline Conserved PrimitiveToConserved(const Primitive &Vc) {
    double kinetic = 0.5 * Vc.rho * (Vc.u*Vc.u + Vc.v*Vc.v);
    double E = Vc.P/(1.4 - 1.0) + kinetic;
    return { Vc.rho, Vc.rho*Vc.u, Vc.rho*Vc.v, E };
}

//-----------------------------------------------------
// Debug Cartesian mesh
//-----------------------------------------------------
inline void convertToCartesianDebug(
    double Ln,
    vector<vector<double>> &xc,
    vector<vector<double>> &yc,
    vector<vector<double>> &Ai,
    vector<vector<double>> &Aj,
    vector<vector<double>> &nxi,
    vector<vector<double>> &nyi,
    vector<vector<double>> &nxj,
    vector<vector<double>> &nyj
) {
    int Ni = imax + 2*ghost;
    int Nj = jmax + 2*ghost;
    double dx = Ln/double(imax);
    double dy = Ln/double(jmax);
    
    xc.assign(Ni, vector<double>(Nj));
    yc.assign(Ni, vector<double>(Nj));
    for(int i=0; i<Ni; i++) {
        for(int j=0; j<Nj; j++) {
            xc[i][j] = (i - ghost + 0.5) * dx;
            yc[i][j] = (j - ghost + 0.5) * dy;
        }
    }
    
    Ai.assign(Ni+1, vector<double>(Nj, dy));
    nxi.assign(Ni+1, vector<double>(Nj, +1.0));
    nyi.assign(Ni+1, vector<double>(Nj,  0.0));
    Aj.assign(Ni,   vector<double>(Nj+1, dx));
    nxj.assign(Ni,   vector<double>(Nj+1,  0.0));
    nyj.assign(Ni,   vector<double>(Nj+1, +1.0));
}

//-----------------------------------------------------
// Cell area via diagonals
//-----------------------------------------------------
inline double computeCellArea(
    int i, int j,
    const vector<vector<double>>& xc,
    const vector<vector<double>>& yc
) {
    double xa = xc[i][j],    ya = yc[i][j];
    double xb = xc[i+1][j],  yb = yc[i+1][j];
    double xc2 = xc[i+1][j+1], yc2 = yc[i+1][j+1];
    double xd = xc[i][j+1],  yd = yc[i][j+1];
    double acx = xc2 - xa,   acy = yc2 - ya;
    double bdx = xb - xd,    bdy = yb - yd;
    return 0.5 * fabs(acx*bdy - acy*bdx);
}

//-----------------------------------------------------
// MMS analytic solution (case 1 only)
//-----------------------------------------------------
struct MmsParams {
    double rho0, rho_x, rho_y, a_rho_x, a_rho_y;
    double u0,   u_x,   u_y,   a_u_x,   a_u_y;
    double v0,   v_x,   v_y,   a_v_x,   a_v_y;
    double p0,   p_x,   p_y,   a_p_x,   a_p_y;
};
static constexpr MmsParams mmsSup = {
    1.0,  0.15, -0.10, 1.0,  0.50,
    800,  50,   -30,  1.5,  0.60,
    800, -75,    40,  0.5,  0.6666667,
    100000, 20000, 50000, 2.0, 1.0
};
static constexpr double PI = 3.14159265358979323846;

inline double rho_mms(double x,double y) {
    double xi=x/L, eta=y/L;
    return mmsSup.rho0 + mmsSup.rho_x*xi + mmsSup.rho_y*eta
         + mmsSup.a_rho_x*sin(PI*xi) + mmsSup.a_rho_y*sin(PI*eta);
}
inline double uvel_mms(double x,double y) {
    double xi=x/L, eta=y/L;
    return mmsSup.u0 + mmsSup.u_x*xi + mmsSup.u_y*eta
         + mmsSup.a_u_x*sin(PI*xi) + mmsSup.a_u_y*sin(PI*eta);
}
inline double vvel_mms(double x,double y) {
    double xi=x/L, eta=y/L;
    return mmsSup.v0 + mmsSup.v_x*xi + mmsSup.v_y*eta
         + mmsSup.a_v_x*sin(PI*xi) + mmsSup.a_v_y*sin(PI*eta);
}
inline double press_mms(double x,double y) {
    double xi=x/L, eta=y/L;
    return mmsSup.p0 + mmsSup.p_x*xi + mmsSup.p_y*eta
         + mmsSup.a_p_x*sin(PI*xi) + mmsSup.a_p_y*sin(PI*eta);
}

//-----------------------------------------------------
// Initialize MMS: fills V and U
//-----------------------------------------------------
void initializeMMS() {
    int Ni = imax + 2*ghost;
    int Nj = jmax + 2*ghost;
    V.assign(Ni, vector<Primitive>(Nj));
    U.assign(Ni, vector<Conserved>(Nj));
    for(int i=0; i<Ni; i++) {
        for(int j=0; j<Nj; j++) {
            double x = x_cell[i][j], y = y_cell[i][j];
            Primitive p{ rho_mms(x,y), uvel_mms(x,y), vvel_mms(x,y), press_mms(x,y) };
            V[i][j] = p;
            U[i][j] = PrimitiveToConserved(p);
        }
    }
    cout << "[DEBUG] MMS init done\n";
}

//-----------------------------------------------------
// Dirichlet BC = exact MMS everywhere in ghosts
//-----------------------------------------------------
void applyBoundaryConditions() {
    int Ni = imax + 2*ghost;
    int Nj = jmax + 2*ghost;
    for(int i=0; i<Ni; i++) {
        for(int j=0; j<Nj; j++) {
            if(i>=ghost && i<ghost+imax && j>=ghost && j<ghost+jmax) continue;
            double x = x_cell[i][j], y = y_cell[i][j];
            Primitive p{ rho_mms(x,y), uvel_mms(x,y), vvel_mms(x,y), press_mms(x,y) };
            V[i][j] = p;
            U[i][j] = PrimitiveToConserved(p);
        }
    }
    cout << "[DEBUG] BC applied\n";
}

//-----------------------------------------------------
// Output initial solution in Tecplot POINT format
//-----------------------------------------------------
void OutputSolution(const string &filename, int iter, bool writeHeader) {
    ofstream file(filename, ios::app);
    if(!file) { cerr << "[ERROR] Cannot open " << filename << "\n"; return; }
    if(writeHeader) {
        file << "TITLE=\"Debug MMS Solution\"\n";
        file << "VARIABLES=\"x\" \"y\" \"rho\" \"u\" \"P\"\n";
    }
    file << "ZONE T=\"Iter " << iter << "\", I=" << imax << ", J=" << jmax << ", F=POINT\n";
    for(int i=ghost; i<ghost+imax; i++) {
        for(int j=ghost; j<ghost+jmax; j++) {
            file << fixed << setprecision(6)
                 << x_cell[i][j] << " ``" << y_cell[i][j] << " "
                 << V[i][j].rho << " " << V[i][j].u << " " << V[i][j].P << "\n";
        }
    }
    cout << "[DEBUG] Wrote solution iter=" << iter << "\n";
}

//-----------------------------------------------------
// main() – debug skeleton
//-----------------------------------------------------
int main() {
    // 1) setup folder
    string outFolder = "OutputFiles";
    if(!fs::exists(outFolder)) {
        fs::create_directory(outFolder);
        cout << "[DEBUG] Created folder " << outFolder << "\n";
    }

    // 2) read header + debug mesh
    const string meshFile = "curv2d9.grd";
    ifstream in(meshFile);
    int nz, kmax;
    in >> nz >> imax >> jmax >> kmax;
    assert(nz == 1 && kmax == 2);
    convertToCartesianDebug(L,
        x_cell, y_cell,
        A_face_i, A_face_j,
        nx_face_i, ny_face_i,
        nx_face_j, ny_face_j);
    cout << "[DEBUG] Debug mesh imax=" << imax << ", jmax=" << jmax << "\n";

    // 3) cell volumes
    vector<vector<double>> cellVol(imax+2*ghost, vector<double>(jmax+2*ghost));
    for(int i=ghost; i<ghost+imax; i++) {
        for(int j=ghost; j<ghost+jmax; j++) {
            cellVol[i][j] = computeCellArea(i, j, x_cell, y_cell);
        }
    }
    cout << "[DEBUG] Cell volumes computed\n";

    // 4) init MMS + BCs
    initializeMMS();
    applyBoundaryConditions();

    // 5) output initial solution
    string solFile = outFolder + "/main_debug_solution.dat";
    OutputSolution(solFile, 0, true);

    // sample check
    cout << "Sample rho at (ghost,ghost)=" << V[ghost][ghost].rho << "\n";
    return 0;
}
