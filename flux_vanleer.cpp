#include "flux_vanleer.hpp"
#include <cmath>

// -- Split Mach number coefficients
inline double Cplus(double M)   { return (M <= -1.0) ? 0.0 : (M < 1.0) ? 0.25 * (M + 1) * (M + 1) : M; }
inline double Cminus(double M)  { return (M >=  1.0) ? 0.0 : (M > -1.0) ? -0.25 * (M - 1) * (M - 1) : M; }
inline double Dplus(double M)   { return (M <= -1.0) ? 0.0 : (M < 1.0) ? Cplus(M) * (-M + 2.0) : 1.0; }
inline double Dminus(double M)  { return (M >=  1.0) ? 0.0 : (M > -1.0) ? Cminus(M) * (-M - 2.0) : 1.0; }

// -- Dot product of velocity with unit normal
inline double dot(const Primitive& V, double nx, double ny) {
    return V.u * nx + V.v * ny;
}

Conserved ComputeVanLeerFlux(
    const Primitive& VL,
    const Primitive& VR,
    const FaceGeometry& face
) {
    const double nx = face.nx;
    const double ny = face.ny;
    const double A  = face.A;

    // LEFT state
    double aL = std::sqrt(gamma_ * VL.P / VL.rho);
    double UL_hat = dot(VL, nx, ny);        // Normal velocity
    double ML = UL_hat / aL;
    double htL = (gamma_ / (gamma_ - 1.0)) * (VL.P / VL.rho) + 0.5 * (VL.u * VL.u + VL.v * VL.v);

    double CpL = Cplus(ML), DpL = Dplus(ML);

    // RIGHT state
    double aR = std::sqrt(gamma_ * VR.P / VR.rho);
    double UR_hat = dot(VR, nx, ny);
    double MR = UR_hat / aR;
    double htR = (gamma_ / (gamma_ - 1.0)) * (VR.P / VR.rho) + 0.5 * (VR.u * VR.u + VR.v * VR.v);

    double CmR = Cminus(MR), DmR = Dminus(MR);

    // -- Convective fluxes
    Conserved F_conv_L, F_conv_R;

    F_conv_L.rho   = VL.rho * aL * CpL;
    F_conv_L.rhou  = VL.rho * aL * CpL * VL.u;
    F_conv_L.rhov  = VL.rho * aL * CpL * VL.v;
    F_conv_L.E     = VL.rho * aL * CpL * htL;

    F_conv_R.rho   = VR.rho * aR * CmR;
    F_conv_R.rhou  = VR.rho * aR * CmR * VR.u;
    F_conv_R.rhov  = VR.rho * aR * CmR * VR.v;
    F_conv_R.E     = VR.rho * aR * CmR * htR;

    // -- Pressure fluxes
    Conserved F_pres_L = {0.0, DpL * VL.P * nx, DpL * VL.P * ny, 0.0};
    Conserved F_pres_R = {0.0, DmR * VR.P * nx, DmR * VR.P * ny, 0.0};

    // -- Total flux
    Conserved flux;
    flux.rho   = (F_conv_L.rho  + F_conv_R.rho)  * A;
    flux.rhou  = (F_conv_L.rhou + F_conv_R.rhou + F_pres_L.rhou + F_pres_R.rhou) * A;
    flux.rhov  = (F_conv_L.rhov + F_conv_R.rhov + F_pres_L.rhov + F_pres_R.rhov) * A;
    flux.E     = (F_conv_L.E    + F_conv_R.E) * A;

    return flux;
}
