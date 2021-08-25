#include <iostream>
#include <iomanip>
#include <tuple>
#include <cmath>
#include <time.h>
#include "amr.hpp"
#include "vis.hpp"
#include "vis_node.hpp"
#include "iflow.hpp"
#include "variable.hpp"
#include "variablenode.hpp"

/** AMGS */
extern "C" {
  #include "/home/numazawa/Solver/AMGS/amgs-s-1.1/src/amgs_c_interface.h"
}
/** AGMG */
extern "C" void dagmg_(int*,double*,int*,int*,double*,double*,int*,int*,int*,int*,double*);

namespace smac
{

  /** Boundary plane */
  void bndpp( int noe ){

    for( int ie = 0; ie < noe; ie++ ){
      for( int i = 0; i < dimn; i++ ){
        if( ele[ie]->ibndn[i] == true ){
          if( ele[ie]->outflown[i] == true ){
            ele[ie]->ppf[i] = 0e0;
            ele[ie]->pf[i]  = 0e0;
          }
          else{
            ele[ie]->ppf[i] = ele[ie]->pp;
            ele[ie]->pf[i]  = ele[ie]->p;
          }
        }
      }
    }
  }
  
  /** Conservation of Momentum Equation*/
  void euler_explicit( double h, int noe ){
  
    for( int ie = 0; ie < noe; ie++ ){
      for( int i = 0; i < dim; i++ ){
        /** Initialize */
        double ct = 0e0, dt = 0e0, pt = 0e0 ;
        /** Pressure term*/
        pt -= ele[ie]->dPdx[i] * ele[ie]->vol / ele[ie]->deno;
        /** Convection term and Diffusion term */
        for( int j = 0; j < dim; j++ ){
          /** Convection term */
          double cttmp1 = 0e0, cttmp2 = 0e0;

          /* At face of East, North and Top*/
          /** 1st-order upwind */
          double tmp1  = ele[ie]->deno * ele[ie]->uf[j][2*j+1];
          if( ele[ie]->ibndn[2*j+1] == false ){
            tmp1 *= std::max( cttmp1, 0e0 ) * ele[ie]->uo[i]
                  + std::min( cttmp1, 0e0 ) * ele[ie]->cnx[2*j+1]->uo[i];
          }
          else{
            tmp1 *= std::max( cttmp1, 0e0 ) * ele[ie]->uo[i]
                  + std::min( cttmp1, 0e0 ) * ele[ie]->uf[i][2*j+1];
          }
          cttmp1 += tmp1 * ele[ie]->sa[2*j+1];
          /** 2nd-order central*/
          cttmp2 += ele[ie]->deno * ele[ie]->uf[i][2*j+1] * ele[ie]->uf[j][2*j+1] * ele[ie]->sa[2*j+1];

          /* At face of West, South and Bottom */
          /** 1st-order upwind */
          tmp1  = ele[ie]->deno * ele[ie]->uf[j][2*j];
          if( ele[ie]->ibndn[2*j] == false ){
            tmp1 *= std::max( cttmp1, 0e0 ) * ele[ie]->cnx[2*j]->uo[i]
                  + std::min( cttmp1, 0e0 ) * ele[ie]->uo[i];
          }
          else{
            tmp1 *= std::max( cttmp1, 0e0 ) * ele[ie]->uf[i][2*j]
                  + std::min( cttmp1, 0e0 ) * ele[ie]->uo[i];
          }
          cttmp1 -= tmp1 * ele[ie]->sa[2*j];
          /** 2nd-order central */
          cttmp2 -= ele[ie]->deno * ele[ie]->uf[i][2*j] * ele[ie]->uf[j][2*j] * ele[ie]->sa[2*j];

          cttmp1 /= ele[ie]->deno;
          cttmp2 /= ele[ie]->deno;
          ct += cttmp1 * 0.0e0 + cttmp2 * 1.0e0;

          /** diffusion term */
          double dttmp = 0e0;

          /* At face of East, North and Top*/
          /* 1st term */
          double tmp = ele[ie]->visc * ele[ie]->dufdx[i][j][2*j+1] * ele[ie]->sa[2*j+1];
          /* 2nd term */
          if ( i == j ) tmp *= 2e0;
          else          tmp += ele[ie]->visc * ele[ie]->dufdx[j][i][2*j+1] * ele[ie]->sa[2*j+1];
          dttmp += tmp;

          /* At face of West, South and Bottom*/
          /* 1st term */
          tmp = ele[ie]->visc * ele[ie]->dufdx[i][j][2*j] * ele[ie]->sa[2*j];
          /* 2nd term */
          if ( i == j ) tmp *= 2e0;
          else          tmp += ele[ie]->visc * ele[ie]->dufdx[j][i][2*j] * ele[ie]->sa[2*j];
          dttmp -= tmp;

          dttmp /= ele[ie]->deno;
          dt += dttmp;
        }

        ele[ie]->u[i] = ele[ie]->uo[i] + ( pt - ct + dt ) / ele[ie]->vol * h;

      }
    }
  }
  
  /** Ordinary Differential Equation of Conservation and Momentum Equation */
  void ode( int noe, bool debug ){
    int status;
    double minu, minv, minw, maxu, maxv, maxw;
    if ( debug ){
      std::cout << std::scientific << std::showpos << std::setprecision(4) \
      << "(*smac::ode) TimeStep: " << delt << " s" << std::endl;
    }
    /** 1st Euler-explicite method*/
    smac::euler_explicit(delt,noe);
    /** 4th Runge-Kutta method*/
    //rk4(delt,noe);
    if ( debug ){
      std::tie(minu,minv,minw) = iflow::minuvw(noe);
      std::tie(maxu,maxv,maxw) = iflow::maxuvw(noe);
      std::cout << std::scientific << std::showpos << std::setprecision(4) \
      << "(*smac::ode) MinVel.: u " << minu << " m/s, v " << minv << " m/s, w " << minw << " m/s" << std::endl;
      std::cout << std::scientific << std::showpos << std::setprecision(4) \
      << "(*smac::ode) MaxVel.: u " << maxu << " m/s, v " << maxv << " m/s, w " << maxw << " m/s" << std::endl;
    }
  }
  
  void calcuf( int noe ){

    for ( int ie = 0; ie < noe; ie++ ){
      for( int i = 0; i < dim; i++ ){
        for( int j = 0; j < dimn; j++ ){
          if( ele[ie]->ibndn[j] == false ){
            double wf = ele[ie]->cnx[j]->xl[0] / ( ele[ie]->xl[0] + ele[ie]->cnx[j]->xl[0] );
            double u1 = ele[ie]->u[i]          + delt / ele[ie]->den         * ele[ie]->dPdx[i]        ;
            double u2 = ele[ie]->cnx[j]->u[i]  + delt / ele[ie]->cnx[j]->den * ele[ie]->cnx[j]->dPdx[i];
            double tmp;
            if( j == 2*i || j == 2*i+1 ){
              tmp  = ( ele[ie]->cnx[j]->p  - ele[ie]->p  ) / ele[ie]->dis_cnx[j];
              if( j == 0 || j == 2 || j == 4 ) tmp  *= -1e0;
            }
            else{
              tmp  = wf * ele[ie]->dPdx[i]  + ( 1e0 - wf ) * ele[ie]->cnx[j]->dPdx[i] ;
            }
            ele[ie]->uf[i][j] = wf * u1  + ( 1e0 - wf ) * u2 - delt / ele[ie]->deno * tmp; 
          }
        }
      }

    }
  }
  
  void coef( int noe ){
  
    for( int ie = 0; ie < noe; ie++ ){
      double nonort_comp = 0e0;
      ele[ie]->sm = 0e0;
      for( int i=0; i < dimn; i++ ){
        if( ele[ie]->ibndn[i] == false ){
          ele[ie]->an[i] = 1e0 / ele[ie]->dis_cnx[i] * ele[ie]->sa[i];
        }
        else{
          if( ele[ie]->outflown[i] == true ) ele[ie]->an[i]  = 2e0 / ele[ie]->xl[0] * ele[ie]->sa[i];
        }
      }
    }

  }
  
  int fvm2crs( int noe, bool AGMG ){
  
    int n = 0;
  
    for( int ie = 0; ie < noe; ie++ ) ele[ie]->pp = 0e0;
  
    if ( AGMG ) crs_index[0] = 1;
    else        crs_index[0] = 0;

    for( int ie = 0; ie < noe; ie++ ){
      /*create coef*/
      x[ie] = 0e0;
      b[ie] = 0e0;
      ele[ie]->ap = 0e0;
      for( int i=0; i < dimn; i++ ){
        if( ele[ie]->ibndn[i] == false ){
          ele[ie]->ap += ele[ie]->an[i];
        }
        else{
          if( ele[ie]->outflown[i] == true ){
            ele[ie]->ap += ele[ie]->an[i];
          }
        }
      }
      /*create crs*/
      if ( AGMG ) crs_col[n] = ele[ie]->key;
      else        crs_col[n] = ele[ie]->key-1;
      crs_val[n] = ele[ie]->ap;
      n++;
      for( int i=0; i < dimn; i++ ){
        if( ele[ie]->ibndn[i] == false ){
          if ( AGMG ) crs_col[n] = ele[ie]->cnx[i]->key;
          else        crs_col[n] = ele[ie]->cnx[i]->key-1;
          crs_val[n] = -ele[ie]->an[i];
          n++;
        }
      }
      if ( AGMG ) crs_index[ie+1] = n+1;
      else        crs_index[ie+1] = n;
    }

    return n;
  
  }
  
  void update_bx( int noe ){
  
    for( int ie = 0; ie < noe; ie++ ){
      b[ie] = ele[ie]->sm;
      x[ie] = 0e0;
      for( int i = 0; i < dim; i++ ){
        /* At face of East, North and Top*/
        b[ie] -= ele[ie]->den * ele[ie]->uf[i][2*i+1] / delt * ele[ie]->sa[2*i+1];
        /* At face of West, South and Bottom*/
        b[ie] += ele[ie]->den * ele[ie]->uf[i][2*i]   / delt * ele[ie]->sa[2*i];
      }
    }

  }
  
  void crs2fvm( int noe ){
    for( int ie = 0; ie < noe; ie++ ) ele[ie]->pp = x[ie];
  }
  
  void solve_s( int noe, int step ){
  
    int ijob, iprint = -6, nrest = 10, iter = 5; 
    double tol = 1.e-10;
    bool AGMG = true;
  
    for( int ie = 0; ie < noe; ie++ ) ele[ie]->pp = 0e0;

    if( step == 1 ){
      noc = 0;
      for( int i = 0; i < mx*my*mz; i++ ){ noc = amr::count_nbr(noc,root[i],noe); }
      /* allocate */
      crs_index = new int[noe+1];
      crs_col   = new int[noc];
      crs_val   = new double[noc];
      b = new double[noe];
      x = new double[noe];
      // fvm -> crs
      noc = smac::fvm2crs(noe,AGMG); 
    }

    smac::update_bx(noe);

    // Matrix solver
    if ( AGMG ) {
      if( step == 1 ){
        ijob = 1;
        dagmg_(&noe,crs_val,crs_col,crs_index,b,x,&ijob,&iprint,&nrest,&iter,&tol);
      }
      ijob = 12;
      dagmg_(&noe,crs_val,crs_col,crs_index,b,x,&ijob,&iprint,&nrest,&iter,&tol);
    }
    else {
      if( step == 1 ){
        amgs_c_parameter_set();
        amgs_c_store_problem_matrix(crs_index, crs_col, crs_val, noe, noc);
        amgs_c_terminal_val_set(tol);
      }
      amgs_c_solver(x, b, noe, &iter, false);
    }
    
    // crs -> fvm
    smac::crs2fvm(noe);
  
    // deallocate
    //delete[] crs_index;
    //delete[] crs_col;
    //delete[] crs_val;
    //delete[] b;
    //delete[] x;
  
  }
  
  void update_p( int noe, double ppref ){
    /** P at cell-center */
    for( int ie = 0; ie < noe; ie++ ) ele[ie]->p += ( ele[ie]->pp - ppref );

  }
  
  /** Pressure Correction */
  void pc( int noe ){
  
    for( int ie = 0; ie < noe; ie++ ){
      for( int i = 0; i < dim; i++ ){
  
        /** Correct velocity at cell-center */
        ele[ie]->u[i] -= delt / ele[ie]->den * ele[ie]->dPPdx[i];
  
        /** Correct velocity at cell-face */
        for( int j = 0; j < dimn; j++ ){
          double tmp = 0e0;
          if( ele[ie]->ibndn[j] == false ){
            if( j == 2*i || j == 2*i+1 ){
              tmp = ( ele[ie]->cnx[j]->pp - ele[ie]->pp ) / ele[ie]->dis_cnx[j];
              if( j == 0 || j == 2 || j == 4 ) tmp  *= -1e0;
            }
            else{
              double wf = ele[ie]->cnx[j]->xl[0] / ( ele[ie]->xl[0] + ele[ie]->cnx[j]->xl[0] );
              tmp  = wf * ele[ie]->dPPdx[i]  + ( 1e0 - wf ) * ele[ie]->cnx[j]->dPPdx[i] ;
            }
            /** At normal plane */
            ele[ie]->uf[i][j] -= delt / ele[ie]->den * tmp;
          }
          /** At outflow boundary */
          if( i == 0 && ele[ie]->ibndn[j] == true && ele[ie]->outflown[j] == true){
            tmp = ( 0e0 - ele[ie]->pp ) / ele[ie]->xl[i] * 2e0;
            if( j == 0 || j == 2 || j == 4 ) tmp  *= -1e0;
            ele[ie]->uf[i][j] -= delt / ele[ie]->den * tmp;
          }
        }
      }
    }
  }
  
  void SMAC( int step, int noe, int non, int outstep, double tol ){

    bool debug = true;
    //================> ok <===================/
    /** Initialize & Visualization */
    if( step == 1 ) {
      iflow::init( noe );
      iflow::bndu( noe, true );
      //ensg_node::ele2node(non);
      //ensg_node::ensg( 0, noe, non );
      ensg::ensg( 0, noe );
    }
    /** CPU time */
    clock_t start = clock();
    /** update velocity, density and viscosity */
    iflow::update( noe );
    ensg_node::ele2node( non );
    iflow::bndu( noe, true ); 
    iflow::gradu( noe ); 
    /** Determine time increment based on Courant number */
    iflow::detdt( noe ); 
    /** p and pp */
    smac::bndpp( noe ); // od
    iflow::gradp_and_pf( noe ); 
    /** Calculate "cell-center velocity*" by Conservation and Momentum Equation */
    smac::ode( noe, false );
    /** Calculate "cell-faced velocity*" by Rhie-Chow interpolation */
    smac::calcuf( noe );
    for( int i = 0; i < 20; i++ ){
      /** boundary */
      iflow::bndu( noe, true );
      /** solve Poisson Equation of pressure increment */
      smac::coef( noe );
      smac::solve_s( noe, step );
      smac::bndpp( noe );
      smac::update_p( noe, ele[0]->pp );
      /** p and pp */
      smac::bndpp( noe );
      iflow::gradp_and_pf( noe ); 
      /** Update velocity */
      smac::pc( noe ); 
      iflow::bndu( noe, false );
      /** Check Conservation Momentum Equation */
      double res = iflow::check_cme( noe, 0e0 ) / reso;
  
      /** CPU time and judge  */
      clock_t end = clock();
      /** Maximum and Minimum of velocities  */
      const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
      if( debug && step%outstep == 0 ) {
        double minu, minv, minw, maxu, maxv, maxw;
        std::tie(minu,minv,minw) = iflow::minuvw(noe);
        std::tie(maxu,maxv,maxw) = iflow::maxuvw(noe);
        std::cout << std::scientific << std::showpos << std::setprecision(4) \
        << "(*smac::smac) MinVel.: u " << minu << " m/s, v " << minv << " m/s, w " << minw << " m/s" << std::endl;
        std::cout << std::scientific << std::showpos << std::setprecision(4) \
        << "(*smac::smac) MaxVel.: u " << maxu << " m/s, v " << maxv << " m/s, w " << maxw << " m/s" << std::endl;
      }
      /** Judge */
      if( std::abs(res) < tol ){
        if( step%outstep == 0 ) {
          std::cout << std::scientific << std::showpos << std::setprecision(4) \
          << "(*smac::smac) "<< i <<" converged!, res: " << res << " kg/s," << delt << " s, time: " << time << std::endl;
          return;
        }
      }
      else{
        if( step%outstep == 0 ) {
          std::cout << std::scientific << std::showpos << std::setprecision(4) \
          << "(*smac::smac) "<< i <<" not converged!, res: " << res << " kg/s, " << delt << " s, time: " << time << std::endl;
        }
      }
    }
  }

}
