#pragma once

//Declare all the physical constant values here

constexpr double Rgas      = 287.0;  // [J/(kg·K)]
constexpr double gamma = 1.4;
constexpr double L = 1.0; //domain length
// We also assume ghost is the number of ghost layers on each side (e.g. 2):
static constexpr int ghost = 2;

