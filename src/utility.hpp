#pragma once

const double PI = 3.141592653589793238462643383279502884;

inline double heaviside(double x){
  return 1. ? x > 0. : 0.;
}
