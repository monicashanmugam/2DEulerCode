// Variables.hpp

#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#pragma once

#include "physicalConstants.hpp"

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

// Declare the globals that were defined in main.cpp:
extern int imax;
extern int jmax;
extern std::vector<std::vector<Conserved>> U;
extern std::vector<std::vector<Primitive>> V;


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
    Ucell.rhov = Vcell.rho * Vcell.v;  // new
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

struct FaceGeometry {
    double nx = 0.0, ny = 0.0;  // Unit normal vector
    double A = 0.0;             // Area (or length in 2D)
};

#endif


