// This header provides two routines for reading a 2D structured curvilinear
// grid from a .grd file (readCurviMeshFromFile) and, in debug mode, for
// overwriting it by a uniform Cartesian grid (convertToCartesianDebug).
//
// Both routines fill the same eight arrays (face-lengths and outward normals),
// so that the rest of your solver can be written identically whether you are
// running on a true curvilinear mesh or on a simple Cartesian test mesh.
// 
// 
#ifndef MESH_UTILS_HPP_
#define MESH_UTILS_HPP_

#include <vector>
#include <fstream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <string>
#include <iostream>
#include "physicalConstants.hpp"
#include "geometry.hpp"
#include "variables.hpp"



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

#endif // MESH_UTILS_HPP_