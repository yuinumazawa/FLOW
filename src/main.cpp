#include <iostream> 
#include <iomanip>
#include <cmath> 
#include <fstream> 
#include <tuple> 
#include <ctime> 
#include "variable.hpp"
#include "variablenode.hpp"
#include "vis.hpp"
#include "vis_node.hpp"
#include "amr.hpp"
#include "iflow.hpp"

/** dimension of AMR */ 
/** if 1 dimension, dim = 1, dimn = 2, dimb = 1(001) */
/** if 2 dimension, dim = 2, dimn = 4, dimb = 3(011) */
/** if 3 dimension, dim = 3, dimn = 8, dimb = 7(111) */
int dim = 2, dimn = 4, dimb = 3;
int mx  = 200, my = 80, mz = 1;
int lvllim= 0;
double xtot = 25e0, ytot = 10e0, ztot = 5e-2;
double delt = 1e20;
double maxdivu = 0e0;
double lim_nonort = 1e0;
double reso = 0e0;
double basecircle;
std::FILE *file;

double *crs_val, *x, *b;
int n, *crs_col, *crs_index, noc;
Element *NIL, **root, **ele, *front;
Node *NILNode, **rootNode, *frontnode;

int main(){
  int noe = 0, nob = 0, non = 0, ts = 100000;
  double time = 0e0;
  int outstep = 100;
  bool out = true;
  std::clock_t cputime = 0e0;
  long double ccputime = 0e0;

  /* Initialization */
  std::cout << "init" << std::endl;
  amr::init(mx*my*mz);
  amr::boundary(mx*my*mz);
  /** mesh convert */
  Element* tmp = NIL;
  noe = 0;
  nob = 0;
  for( int i = 0; i < mx*my*mz; i++ ) std::tie(tmp,noe) = amr::zordering(root[i],tmp,noe); 
  for( int i = 0; i < mx*my*mz; i++ ) std::tie(tmp,nob) = amr::zordering_bnd(root[i],tmp,nob); 
  amr::convert_amr2fvm(noe,nob);
  non = ensg_node::zordering_node( noe );

  for ( int step = 1; step <= ts; step++ ){
    cputime = std::clock();
    if( step%outstep == 0 ) {
      out = true;
      std::cout << "(*main) step: " << step << std::endl ;
    }
    else out = false;

    /** Convection */
    iflow::IncompressibleFlow(step, noe, non, outstep, 1e-6);
    cputime = 1000.0 * (std::clock() - cputime) / CLOCKS_PER_SEC;
    ccputime += cputime; 
    time += delt;

    /** Visualization */
    if( step%outstep == 0 ){
      ensg::ensg( step/outstep, noe );
      //ensg_node::ele2node( non );
      //ensg_node::ensg( step/outstep, noe, non );
      std::cout << "(*main) Cumulative time: " << time << " s, delt: " << delt << " s"<< std::endl;
      std::cout << "(*main) Cumulative CPU time: " << ccputime / 1000.0 << " s, CPU time per step: " 
      << cputime / 1000.0 << " s" << std::endl << std::endl ;
      
    }
  }
  amr::convert_fvm2amr(noe,nob);
  return 0;
}
