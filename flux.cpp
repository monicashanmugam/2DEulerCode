#include "flux.hpp"
#include "flux_vanleer.hpp"
#include "flux_roe.hpp"

Conserved ComputeNumericalFlux(
    const Primitive& VL, const Primitive& VR,
    const FaceGeometry& face, FluxScheme scheme
) {
    switch (scheme) {
        case FluxScheme::VanLeer:
            return ComputeVanLeerFlux(VL, VR, face);
        case FluxScheme::Roe:
            return ComputeRoeFlux(VL, VR, face);
        default:
            throw std::runtime_error("Unknown flux scheme.");
    }
}
