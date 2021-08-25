#ifndef SMAC_HPP
#define SMAC_HPP
#include <tuple>
#include "variable.hpp"

namespace smac {
  
  extern void bndpp(int);
  
  extern void euler_explicit(double,int);
  extern void ode(int,bool);
  
  extern void coef(int);
  extern int  fvm2crs(int,bool);
  extern void update_bx(int);
  extern void crs2fvm(int);
  extern void solve_agmg_s(int,bool);
  extern void update_p(int,double);
  extern void pc(int);
  
  extern void SMAC(int,int,int,int,double);
  
}

#endif
