// flux.hpp
#pragma once
#include "variables.hpp"

enum class FluxScheme { VanLeer, Roe };

Conserved ComputeNumericalFlux(
    const Primitive& VL, const Primitive& VR,
    const FaceGeometry& face, FluxScheme scheme
);
