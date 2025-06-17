#include "fluxSchemes.hpp"
#include <cmath>

Conserved ComputeFluxAcrossFace(
    const Primitive& VL,
    const Primitive& VR,
    const FaceGeometry& face,
    double gamma,
    FluxScheme scheme
) {
    switch (scheme) {
        case ROE:
            return RoeFlux(VL, VR, face, gamma);
        case VANLEER:
        default:
            return VanLeerFlux(VL, VR, face, gamma);
    }
}

// Helper: compute normal velocity and tangential splitting
Conserved VanLeerFlux(
    const Primitive& VL,
    const Primitive& VR,
    const FaceGeometry& face,
    double gamma
) {
    // Rotate left and right states into local face-normal frame
    double unL = VL.u * face.nx + VL.v * face.ny;
    double unR = VR.u * face.nx + VR.v * face.ny;

    double H_L = (VL.P / (gamma - 1.0) + 0.5 * VL.rho * (VL.u * VL.u + VL.v * VL.v) + VL.P) / VL.rho;
    double H_R = (VR.P / (gamma - 1.0) + 0.5 * VR.rho * (VR.u * VR.u + VR.v * VR.v) + VR.P) / VR.rho;

    double aL = std::sqrt(gamma * VL.P / VL.rho);
    double aR = std::sqrt(gamma * VR.P / VR.rho);

    // Compute mass flux using Van Leer's flux vector splitting
    auto flux = Conserved{};

    auto f_plus = [&](const Primitive& V, double un, double a, double H) -> Conserved {
        Conserved F;
        if (un >= 0) {
            F.rho  = V.rho * un;
            F.rhou = F.rho * V.u;
            F.rhov = F.rho * V.v;
            F.E    = F.rho * H;
        }
        return F;
    };

    auto f_minus = [&](const Primitive& V, double un, double a, double H) -> Conserved {
        Conserved F;
        if (un < 0) {
            F.rho  = V.rho * un;
            F.rhou = F.rho * V.u;
            F.rhov = F.rho * V.v;
            F.E    = F.rho * H;
        }
        return F;
    };

    Conserved Fplus  = f_plus(VL, unL, aL, H_L);
    Conserved Fminus = f_minus(VR, unR, aR, H_R);

    flux = Fplus + Fminus;

    // Multiply by face area (for integration)
    flux = flux * face.A;
    return flux;
}
