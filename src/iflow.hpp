#ifndef UFLOW_HPP
#define UFLOW_HPP
#include <tuple>
#include "variable.hpp"

namespace iflow
{
  extern double sign(double);
  extern void inflow_bndu(int);
  extern void init(int);
  extern void update(int);
  extern void bndu(int,bool);
  extern void detdt(int);
  extern std::tuple<double,double,double>maxuvw(int);
  extern std::tuple<double,double,double>minuvw(int);
  extern void update(int);
  extern void gradu(int);
  extern void gradp_and_pf(int);
  extern double check_cme(int,double);
  extern void IncompressibleFlow(int,int,int,int,double);
}

#endif
