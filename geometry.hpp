// geometry.hpp
#pragma once

#include <vector>

//-------------------------------------------------------------
// Cell center coordinates (i,j)
//-------------------------------------------------------------
extern std::vector<std::vector<double>> x_cell;  // x-coordinate of cell center
extern std::vector<std::vector<double>> y_cell;  // y-coordinate of cell center

//-------------------------------------------------------------
// Face areas and unit normals (i = vertical faces, j = horizontal faces)
//-------------------------------------------------------------

// Face areas (lengths in 2D)
extern std::vector<std::vector<double>> A_face_i;  // face between (i-1,j) and (i,j)
extern std::vector<std::vector<double>> A_face_j;  // face between (i,j-1) and (i,j)

// Outward unit normals (pointing *into* the cell at (i,j))
extern std::vector<std::vector<double>> nx_face_i; // x-component of i-face normal
extern std::vector<std::vector<double>> ny_face_i; // y-component of i-face normal

extern std::vector<std::vector<double>> nx_face_j; // x-component of j-face normal
extern std::vector<std::vector<double>> ny_face_j; // y-component of j-face normal
