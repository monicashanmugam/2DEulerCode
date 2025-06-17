// fluxSchemes.hpp
#pragma once
#include "variables.hpp"

// Structure to hold face-normal and area
struct FaceGeometry {
    double nx, ny; // Unit normal vector (already in your variables.hpp)
    double A;      // Face area (1.0 in 2D)
};

// Enum for selecting scheme
enum FluxScheme { VANLEER, ROE };

// Master interface: returns flux at a face
Conserved ComputeFluxAcrossFace(
    const Primitive& VL,
    const Primitive& VR,
    const FaceGeometry& face,
    double gamma,
    FluxScheme scheme = VANLEER
);

// Low-level functions
Conserved VanLeerFlux(
    const Primitive& VL,
    const Primitive& VR,
    const FaceGeometry& face,
    double gamma
);

Conserved RoeFlux(
    const Primitive& VL,
    const Primitive& VR,
    const FaceGeometry& face,
    double gamma
);
