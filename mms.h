#ifndef MMS_H
#define MMS_H

#pragma once
#include "your_structs.h"  // make sure this includes Primitive definition

namespace MMS {
  void setCase(int mmsCase_, double L_);
  Primitive evaluatePrimitive(double x, double y);
}


namespace MMS {

  // POD to hold all amplitudes for one MMS case
  struct MmsParams {
    double rho0, rho_x, rho_y, a_rho_x, a_rho_y;
    double u0,   u_x,   u_y,   a_u_x,   a_u_y;
    double v0,   v_x,   v_y,   a_v_x,   a_v_y;
    double p0,   p_x,   p_y,   a_p_x,   a_p_y;
  };

  // Call once after mesh geometry is ready:
  //   mmsCase = 1 → supersonic; 2 → subsonic
  //   L       = domain length (e.g. 1.0)
  void initializeMMS(int mmsCase, double L);

  // The analytic source term for ∂ρ/∂t
  double rmassconv(double x, double y);
  double xMtmConv(double x, double y);
  double yMtmConv(double x, double y);
  double EnergyConv(double x, double y);

} // namespace MMS

#endif // MMS_H
