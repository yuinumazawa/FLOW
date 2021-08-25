#include <iostream>
#include <iomanip>
#include <tuple>
#include <cmath>
#include <time.h>
#include "amr.hpp"
#include "smac.hpp"
#include "variable.hpp"

namespace iflow
{

  /** Boundary */
  void inflow_bndu( int noe ){
    reso = 0e0;
    for( int i = 0; i < noe; i++ ){
      for ( int j = 0; j < 6; j++ ){
        if( ele[i]->ibndn[j] == true ){
          if( j == 0 ){
            ele[i]->uf[0][j] = 1e0; // m/s 
            reso += ele[i]->den * ele[i]->uf[0][j] * ele[i]->xl[0] * ele[i]->xl[2];
          }
        }
      }
    }
  }
  
  void init( int noe ){
    for ( int ie = 0; ie < noe; ie++ ){
      ele[ie]->den  = 1.0e0;
      ele[ie]->deno = 1.0e0;
      ele[ie]->visc = 1.0e-2; // Re = 1/ visc -> 1.0e-2 at Re = 100, 5.0e-2 at Re = 20
      ele[ie]->pp   = 0e0;
      ele[ie]->p    = 0e0;
      for( int id = 0; id < 3; id++ ){
        ele[ie]->u[id]     = 0e0;
        ele[ie]->uo[id]    = 0e0;
        ele[ie]->dPdx[id]  = 0e0;
        ele[ie]->dPPdx[id] = 0e0;
        for( int k = 0; k < 6; k++ ){
          ele[ie]->uf[id][k] = 0e0;
          ele[ie]->pf[k]     = 0e0;
          ele[ie]->ppf[k]    = 0e0;
        }
      }
    }
    iflow::inflow_bndu(noe);
  }

  void bndu( int noe, bool iout ){

    for( int ie = 0; ie < noe; ie++ ){
      for ( int i = 0; i < dim; i++ ){
        for ( int j = 0; j < dimn; j++ ){
          if( ele[ie]->ibndn[j] == true ){

            /** non-slip boundary ( wall ) */
            if ( ele[ie]->walln[j] == true ){ 
              ele[ie]->uf[i][j] = 0e0;
              for ( int l = 0; l < dimn; l++ ){
                if( j == 2*l || j == 2*l+1 ) {
                  ele[ie]->dufdx[i][l][j] = (ele[ie]->u[i]-ele[ie]->uf[i][j])/ele[ie]->xl[l]*2e0;
                  if( j%2 != 0 ) ele[ie]->dufdx[i][l][j] *= -1e0;
                }
                else {
                  ele[ie]->dufdx[i][l][j] = 0e0;
                }
              }
            }
            /** moving wall */
            else if ( ele[ie]->movingwalln[j] == true ){ 
              /** ele[ie]->uf[i][j][k]: set in inflowu */
              if( i == 1 ) ele[ie]->uf[i][j] = 1e0;
              else ele[ie]->uf[i][j]         = 0e0;
              for ( int l = 0; l < 3; l++ ){
                if( j == 2*l || j == 2*l+1 ) {
                  ele[ie]->dufdx[i][l][j] = (ele[ie]->u[i]-ele[ie]->uf[i][j])/ele[ie]->xl[l]*2e0;
                  if( j%2 != 0 ) ele[ie]->dufdx[i][l][j] *= -1e0;
                }
                else {
                  ele[ie]->dufdx[i][l][j] = 0e0;
                }
              }
            }
            /** inflow boundary */
            else if ( ele[ie]->inflown[j] ){  
              if( i == 0 ) {
                ele[ie]->uf[i][j] = 1e0;
                for ( int l = 0; l < 3; l++ ){
                  if( j == 2*l || j == 2*l+1 ) {
                    ele[ie]->dufdx[i][l][j] = (ele[ie]->u[i]-ele[ie]->uf[i][j])/ele[ie]->xl[l]*2e0;
                    if( j%2 != 0 ) ele[ie]->dufdx[i][l][j] *= -1e0;
                  }
                  else ele[ie]->dufdx[i][l][j] = 0e0;
                }
              }
              else {
                ele[ie]->uf[i][j] = 0e0;
                for ( int l = 0; l < 3; l++ ){
                  if( j == 2*l || j == 2*l+1 ) {
                    ele[ie]->dufdx[i][l][j] = (ele[ie]->u[i]-ele[ie]->uf[i][j])/ele[ie]->xl[l]*2e0;
                    if( j%2 != 0 ) ele[ie]->dufdx[i][l][j] *= -1e0;
                  }
                  else ele[ie]->dufdx[i][l][j] = 0e0;
                }
              }
              
            }
            /** outflow boundary */
            else if ( ele[ie]->outflown[j] ){ 
              if( j == 2*i || j == 2*i+1 ) {
                if( iout ){
                  ele[ie]->uf[i][j] = ele[ie]->u[i];
                }
                for ( int l = 0; l < 3; l++ ){
                  if( j == 2*l || j == 2*l+1 ) {
                    ele[ie]->dufdx[i][l][j] = (ele[ie]->u[i]-ele[ie]->uf[i][j])/ele[ie]->xl[l]*2e0;
                    if( j%2 != 0 ) ele[ie]->dufdx[i][l][j] *= -1e0;
                  }
                  else ele[ie]->dufdx[i][l][j] = ele[ie]->dudx[i][l];
                }
              }
              else{
                ele[ie]->uf[i][j] = ele[ie]->u[i];
                ele[ie]->dufdx[i][0][j] = 0e0;
                ele[ie]->dufdx[i][1][j] = 0e0;
                ele[ie]->dufdx[i][2][j] = 0e0;
              }
            }
            /** slip boundary */
            else if ( ele[ie]->slipn[j] ){ 
              if( j == 2*i || j == 2*i+1 ) {
                ele[ie]->uf[i][j] = 0e0;
                for ( int l = 0; l < 3; l++ ){
                  if( j == 2*l || j == 2*l+1 ) {
                    ele[ie]->dufdx[i][l][j] = (ele[ie]->u[i]-ele[ie]->uf[i][j])/ele[ie]->xl[l]*2e0;
                    if( j%2 != 0 ) ele[ie]->dufdx[i][l][j] *= -1e0;
                  }
                  else ele[ie]->dufdx[i][l][j] = 0e0;
                }
              }
              else{
                ele[ie]->uf[i][j] = ele[ie]->u[i];
                for ( int l = 0; l < 3; l++ ){
                  if( j == 2*l || j == 2*l+1 ) {
                    ele[ie]->dufdx[i][l][j] = 0e0; 
                  }
                  else ele[ie]->dufdx[i][l][j] = 0e0;
                }
              }
            }

          }
        }
      }
    }

  }

  void detdt( int noe ){
    // Cr = u * dt / dx
    double cr = 0.4e0;
    for( int ie = 0; ie < noe; ie++ ){
      for( int id = 0; id < dim; id++ ){
        double dtmp = cr * ele[ie]->xl[id] / std::max( std::abs( ele[ie]->uo[id]) , 1e-20  );
        delt = std::min( delt, dtmp );
        for ( int j = 0; j < 6; j++ ){
          double dtmp = cr * ele[ie]->xl[id] / std::max( std::abs( ele[ie]->uf[id][j]) , 1e-20  );
          delt = std::min( delt, dtmp );
        }
      }
    }
  }
  
  std::tuple<double,double,double> maxuvw( int noe ){
    double maxu[3];
    maxu[0]=0e0, maxu[1]=0e0, maxu[2]=0e0;
    for( int ie = 0; ie < noe; ie++ ){
      for( int j = 0; j < 3; j++ ){
        maxu[j] = std::max( ele[ie]->u[j], maxu[j] );
      }
    }
    return std::forward_as_tuple(maxu[0],maxu[1],maxu[2]);
  }
  
  std::tuple<double,double,double> minuvw( int noe ){
    double minu[3];
    minu[0]=1e20, minu[1]=1e20, minu[2]=1e20;
    for( int ie = 0; ie < noe; ie++ ){
      for( int j = 0; j < 3; j++ ){
        minu[j] = std::min( ele[ie]->u[j], minu[j] );
      }
    }
    return std::forward_as_tuple(minu[0],minu[1],minu[2]);
  }

  /** update */
  void update( int noe ){
    for ( int ie = 0; ie < noe; ie++ ){
      ele[ie]->deno = ele[ie]->den;
      for ( int i = 0; i < dim; i++ ){
        ele[ie]->uo[i] = ele[ie]->u[i];
      }
    }
  }
  
  /** Gradient of u */ 
  void gradu( int noe ){
  
    /** Gradient of u at cell center*/
    for( int ie = 0; ie < noe; ie++ ){
      for( int i = 0; i < dim; i++ ){
        for( int j = 0; j < dim; j++ ){
          ele[ie]->dudx[i][j]  = 0e0;
          if( i == j ) continue;
          else{
            /* At face of East, North and Top*/
            ele[ie]->dudx[i][j]  += ele[ie]->uf[i][2*j+1] * ele[ie]->sa[2*j+1];
            /* At face of West, South and Bottom*/
            ele[ie]->dudx[i][j]  -= ele[ie]->uf[i][2*j]   * ele[ie]->sa[2*j];
          }
          ele[ie]->dudx[i][j] /= ele[ie]->vol;
        }
      }
    }

    /** Gradient of u at cell face */
    for( int ie = 0; ie < noe; ie++ ){
      for( int i = 0; i < dim; i++ ){
        for( int j = 0; j < dim; j++ ){
          for( int id = 0; id < dimn; id++ ){
            if( ele[ie]->ibndn[id] == false ){
              if( ele[ie]->cnx[id]->lvl == ele[ie]->lvl  ){
                if( id == 2*j+1 || id == 2*j ){
                  ele[ie]->dufdx[i][j][id] = ( ele[ie]->cnx[id]->uo[i] - ele[ie]->uo[i] ) / ele[ie]->dis_cnx[id]; 
                  if( id == 2*j ) ele[ie]->dufdx[i][j][id] *= -1e0;
                }
                else{
                  ele[ie]->dufdx[i][j][id] = ( ( ele[ie]->cnx[id]->uf[i][2*j+1] + ele[ie]->uf[i][2*j+1] ) / 2e0 
                                             - ( ele[ie]->cnx[id]->uf[i][2*j]   + ele[ie]->uf[i][2*j]   ) / 2e0 ) 
                                               / ele[ie]->xl[j]; 
                }
              }
            }
          }
        }
      }
    }
  
  }

  /** Gradient of p and pp */ 
  void gradp_and_pf( int noe ){
    /** P at cell-face of normal boundary*/
    for( int ie = 0; ie < noe; ie++ ){
      for( int j = 0; j < dimn; j++ ){
        if( ele[ie]->ibndn[j] == false ){
          double wf = ele[ie]->cnx[j]->xl[0] / ( ele[ie]->xl[0] + ele[ie]->cnx[j]->xl[0] );
          ele[ie]->pf[j]  =  wf * ele[ie]->p  + ( 1e0 - wf ) * ele[ie]->cnx[j]->p ;
          ele[ie]->ppf[j] =  wf * ele[ie]->pp + ( 1e0 - wf ) * ele[ie]->cnx[j]->pp;
        }
      }
    }

    /** Gradient of P and PP at cell center*/
    for( int ie = 0; ie < noe; ie++ ){
      for( int i = 0; i < dim; i++ ){
        ele[ie]->dPdx[i]  = 0e0;
        ele[ie]->dPPdx[i] = 0e0;
        /* At face of East, North and Top*/
        ele[ie]->dPdx[i]  += ele[ie]->pf[2*i+1]  * ele[ie]->sa[2*i+1];
        ele[ie]->dPPdx[i] += ele[ie]->ppf[2*i+1] * ele[ie]->sa[2*i+1];
        /* At face of West, South and Bottom*/
        ele[ie]->dPdx[i]  -= ele[ie]->pf[2*i]  * ele[ie]->sa[2*i];
        ele[ie]->dPPdx[i] -= ele[ie]->ppf[2*i] * ele[ie]->sa[2*i];

        ele[ie]->dPdx[i]  /= ele[ie]->vol;
        ele[ie]->dPPdx[i] /= ele[ie]->vol;
      }
    }
  
  }
  
  double check_cme( int noe, double res ){
    for( int ie = 0; ie < noe; ie++ ){
      double restmp = 0e0;
      for( int i = 0; i < 3; i++ ){
        /* At face of East, North and Top*/
        restmp += ele[ie]->den * ele[ie]->uf[i][2*i+1] * ele[ie]->sa[2*i+1];
        /* At face of West, South and Bottom*/
        restmp -= ele[ie]->den * ele[ie]->uf[i][2*i]   * ele[ie]->sa[2*i];
      }
      restmp /= ele[ie]->vol;
      ele[ie]->res = std::abs( restmp*ele[ie]->vol ) ;
      res += std::abs( ele[ie]->res );
    }
    return res;
  }

  void IncompressibleFlow( int step, int noe, int non, int outstep, double tol ){
    smac::SMAC(step,noe,non,outstep,tol);
  }

}
