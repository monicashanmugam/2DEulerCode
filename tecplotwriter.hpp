#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <variables.hpp>


// Declare the globals that were defined in main.cpp:
extern int imax;
extern int jmax;
extern std::vector<std::vector<Conserved>> U;
extern std::vector<std::vector<Primitive>> V;

// ----------------------------------------------------------------------------
// Write the Tecplot “VARIABLES=…” header (call once at the very beginning)
// ----------------------------------------------------------------------------
void writeTecplotHeader(std::ofstream &out) {
    out << std::fixed << std::setprecision(6);
    out << "variables=\"x(m)\" \"y(m)\" \"rho(kg/m^3)\" \"u(m/s)\" "
           "\"v(m/s)\" \"press(N/m^2)\"" 
        << "\n";
}

// ----------------------------------------------------------------------------
// Write a single “zone” block (e.g. after each 100 iterations, or at n=0)
//   n         = zone index (an integer, will show up as T="<n>")
//   x_cell, y_cell  = 2D arrays of size [imax+2*ghost][jmax+2*ghost] (cell centers)
//   V             = 2D array of Primitive, size [imax+2*ghost][jmax+2*ghost]
// ----------------------------------------------------------------------------
void writeTecplotZone(
    std::ofstream &out,
    int zoneIndex,
    const std::vector<std::vector<double>> &x_cell,
    const std::vector<std::vector<double>> &y_cell,
    const std::vector<std::vector<Primitive>> &V
) {
    // 1) Zone header
    out << "zone T=\"" << zoneIndex << "\" ";
    out << "I=" << imax << " J=" << jmax << "\n";
    out << "DATAPACKING=POINT\n";

    // 2) Loop in “j, i” order exactly as in the Fortran example
    for (int j = 1; j <= jmax; ++j) {
      for (int i = 1; i <= imax; ++i) {
        // Remember: in your arrays, you have ghost layers at 0..ghost-1
        // so the “first physical cell” is at index (ghost, ghost).
        int gi = i - 1 + ghost;
        int gj = j - 1 + ghost;
        const Primitive &p = V[gi][gj];
        out << x_cell[gi][gj] << " "
            << y_cell[gi][gj] << " "
            << p.rho << " "
            << p.u   << " "
            << p.v   << " "
            << p.P
            << "\n";
      }
    }
}