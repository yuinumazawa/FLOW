#ifndef AMR_HPP
#define AMR_HPP
#include <tuple> 
#include "variable.hpp"
#include "variablenode.hpp"

namespace amr
{
  extern int convert_amr2fvm(int,int);
  extern void convert_fvm2amr(int,int);
  extern int numbering(int, Element*);
  extern int count_nbr(int, Element*,int);
  /** Initialization */
  extern void init(int);
  extern void boundary(int);
  /** Filling curve */
  extern void cr8_cnx(Element*);
  extern std::tuple<Element*,int> zordering(Element*,Element*,int);
  extern std::tuple<Element*,int> zordering_bnd(Element*,Element*,int);
}
#endif
