#pragma once
#include "varialbes.hpp"

Conserved ComputeVanLeerFlux(
    const Primitive& VL,
    const Primitive& VR,
    const FaceGeometry& face
);
