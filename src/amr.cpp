#include <cmath>
#include <iostream>
#include <tuple> 
#include "variable.hpp"
#include "variablenode.hpp"

double costheta = (3e0/sqrt(1e1));

namespace amr
{

  void convert_amr2fvm(int noe, int nob ){
  
    Element* tmp = front;
    ele = new Element*[noe+nob];
    tmp = front;
    for( int ie = 0; ie < noe+nob; ie++ ){
      ele[ie] = tmp;
      tmp = tmp->next;
    }
    for( int ie = 0; ie < noe+nob; ie++ ){
      for( int i = 0; i < dimn; i++ ){
        if( !ele[ie]->ibndn[i] ){
          ele[ie]->cnx[i] = ele[ele[ie]->cnx[i]->key-1];
        }
      }
    }
  }
  
  void convert_fvm2amr( int noe , int nob ){
    Element* tmp = front;
    for( int ie = 0; ie < noe+nob; ie++ ){
      tmp = ele[ie];
      tmp = tmp->next;
    }
    delete [] ele;
  }
 

  void boundary( int noe ){

    for( int ie = 0; ie < noe; ie++ ){

      for( int i = 0; i < dimn; i++ ){
        if( root[ie]->ibndn[i] == true ) {
          if( root[ie]->outside[i] == true ) {
            if( i == 0 ) root[ie]->inflown[i] = true;
            else if( i == 1 ) root[ie]->outflown[i] = true;
            else if( i == 2 || i == 3 ) root[ie]->slipn[i] = true;
            //else if( i == 2 || i == 3 ) root[ie]->walln[i] = true;
            //if( i == 0 ) root[ie]->movingwalln[i] = true;
          }
        }
      }
      
      // karman
      double r  = std::pow(root[ie]->xc[0]-10e0,2e0)+std::pow(root[ie]->xc[1],2e0);
      if ( r < 1.0 )  {
        root[ie]->ibnd = true;
        root[ie]->wall = true;
        if ( (dimb&1) == 1 ){
          root[ie]->nbr[0]->ibndn[1] = true;
          root[ie]->nbr[0]->walln[1] = true;
          root[ie]->nbr[1]->ibndn[0] = true;
          root[ie]->nbr[1]->walln[0] = true;
        }
        if ( (dimb&2) == 2 ){
          root[ie]->nbr[2]->ibndn[3] = true;
          root[ie]->nbr[2]->walln[3] = true;
          root[ie]->nbr[3]->ibndn[2] = true;
          root[ie]->nbr[3]->walln[2] = true;
        }
        if ( (dimb&4) == 4 ){
          root[ie]->nbr[4]->ibndn[5] = true;
          root[ie]->nbr[4]->walln[5] = true;
          root[ie]->nbr[5]->ibndn[4] = true;
          root[ie]->nbr[5]->walln[4] = true;
        }
      }

    }

  }

  int numbering( int n, Element* node ){

    if(node->leaf == true ){
      n++;
      node->key = n;
    }

    return n;
  }
  
  int count_nbr( int n, Element* node, int noe ){

    if( node->leaf == true ){
      n++;
      for ( int i = 0; i < 2*dim; i++ ){
        if( node->ibndn[i] == false ) n++;
      }
    }

    return n;
  }
  
  void cr8_cnx( Element* node ){
    // neiboring
    for( int i = 0; i < 6; i++ ){
      node->cnx[i] = NIL;
      node->dis_cnx[i] = 0e0;
      node->num_cnx[i] = 0;
    }
    for( int i = 0; i < dimn; i++ ){
      if( node->nbr[i] != NIL ){ 
        node->cnx[i]  = node->nbr[i];
        if( node->cnx[i]->ibnd == true ) node->ibndn[i] = true;
        if( node->cnx[i]->wall == true ) node->walln[i] = true;
        if( i == 0 || i == 1 ){
          node->dis_cnx[i] = node->xl[0];
          node->sa[i] = node->xl[1] * node->xl[2];
        }
        else if( i == 2 || i == 3 ){
          node->dis_cnx[i] = node->xl[1];
          node->sa[i] = node->xl[0] * node->xl[2];
        }
        else if( i == 4 || i == 5 ){
          node->dis_cnx[i] = node->xl[2];
          node->sa[i] = node->xl[0] * node->xl[1];
        }
      }
    }
  
  }

  void cr8_bnd( int id, Element* node, bool child, bool outside  ){
    Element *tmpnode;
    tmpnode = new Element;
    node->nbr[id] = tmpnode;
    node->ibndn[id] = true;
    
    tmpnode->ibnd = true;

    /** base */
    if ( !child && outside ) {
      node->outside[id] = true;
    }
    /** child */
    else {
      if ( child && node->parent->outside[id] == true ) {
        node->outside[id] = true;
      }
    }
  }
  
  // initialize
  int init(int key) {
    root  = new Element*[mx*my*mz];
    int n = 0;
    for ( int k = 0; k < mz; k++ ){
      for ( int j = 0; j < my; j++ ){
        for ( int i = 0; i < mx; i++ ){
          /** */
          root[n] = new Element;
          /** Initialize */
          root[n]->key = key;
          /** Geometry */
          root[n]->xi[0] = i;
          root[n]->xi[1] = j;
          root[n]->xi[2] = k;
          root[n]->xl[0] = xtot/double(mx);
          root[n]->xl[1] = ytot/double(my);
          root[n]->xl[2] = ztot/double(mz);
          root[n]->xc[0] = root[n]->xl[0]/2e0+double(i)*root[n]->xl[0];//-double(5e-1)*xtot;
          root[n]->xc[1] = root[n]->xl[1]/2e0+double(j)*root[n]->xl[1]-double(5e-1)*ytot;
          root[n]->xc[2] = root[n]->xl[2]/2e0+double(k)*root[n]->xl[2]-double(5e-1)*ztot;
          for( int h = 0; h < 8; h++ ){
            if ( (h&1) == 1 ) root[n]->xg[0][h] = root[n]->xc[0] + root[n]->xl[0]/2e0;
            else              root[n]->xg[0][h] = root[n]->xc[0] - root[n]->xl[0]/2e0;
            if ( (h&2) == 2 ) root[n]->xg[1][h] = root[n]->xc[1] + root[n]->xl[1]/2e0;
            else              root[n]->xg[1][h] = root[n]->xc[1] - root[n]->xl[1]/2e0;
            if ( (h&4) == 4 ) root[n]->xg[2][h] = root[n]->xc[2] + root[n]->xl[2]/2e0;
            else              root[n]->xg[2][h] = root[n]->xc[2] - root[n]->xl[2]/2e0;
            root[n]->nbrnode[h] = NILNode;
          }
          root[n]->vol = root[n]->xl[0]*root[n]->xl[1]*root[n]->xl[2];
          /** Offspring relationship */
          root[n]->parent = root[n];  
          for ( int h = 0; h < dimn; h++ ) root[n]->child[h] = NIL;
          for ( int h = 0; h < 6; h++ )    root[n]->nbr[h] = NIL;
          for ( int h = 0; h < 6; h++ ) {
            root[n]->cnx[h] = NIL;
            root[n]->dis_cnx[h] = 0e0;
            root[n]->ibndn[h] = false;
            root[n]->slipn[h] = false;
            root[n]->walln[h] = false;
            root[n]->movingwalln[h] = false;
            root[n]->outside[h] = false;
            root[n]->inflown[h] = false;
            root[n]->outflown[h] = false;
          }
          root[n]->lvl = 0;
          root[n]->prev = NIL;
          root[n]->next = NIL;
          /** */
          root[n]->val  =  0e0;
          root[n]->valo =  0e0;
          root[n]->u[0]  = 0e0;
          root[n]->u[1]  = 0e0;
          root[n]->u[2]  = 0e0;
          root[n]->uo[0] = 0e0;
          root[n]->uo[1] = 0e0;
          root[n]->uo[2] = 0e0;
          root[n]->den   = 0e0;
          root[n]->deno  = 0e0;
          root[n]->p     = 0e0;
          root[n]->pp    = 0e0;
          root[n]->gamm  = 0e0;
          root[n]->visc  = 0e0;
          root[n]->vist  = 0e0;
          /** */
          root[n]->ap  = 0e0;
          root[n]->apo = 0e0;
          /** */
          n++;
          key++;
        }
      }
    }
    /** Neighboring cell */
    for( int i=0; i < mx*my*mz; i++ ){
      if( root[i]->xi[2] != 0    && (dimb&4) == 4 ) root[i]->nbr[4] = root[i-mx*my]; 
      else                                          amr::cr8_bnd(4,root[i],false,true);
      if( root[i]->xi[1] != 0    && (dimb&2) == 2 ) root[i]->nbr[2] = root[i-mx]; 
      else                                          amr::cr8_bnd(2,root[i],false,true);
      if( root[i]->xi[0] != 0    && (dimb&1) == 1 ) root[i]->nbr[0] = root[i-1]; 
      else                                          amr::cr8_bnd(0,root[i],false,true);
      if( root[i]->xi[0] != mx-1 && (dimb&1) == 1 ) root[i]->nbr[1] = root[i+1]; 
      else                                          amr::cr8_bnd(1,root[i],false,true);
      if( root[i]->xi[1] != my-1 && (dimb&2) == 2 ) root[i]->nbr[3] = root[i+mx]; 
      else                                          amr::cr8_bnd(3,root[i],false,true);
      if( root[i]->xi[2] != mz-1 && (dimb&4) == 4 ) root[i]->nbr[5] = root[i+mx*my]; 
      else                                          amr::cr8_bnd(5,root[i],false,true);
    }

    /** Neighboring cell from Node */
    rootNode = new Node*[(mx+1)*(my+1)*(mz+1)];
    n = 0;
    for ( int k = 0; k < mz+1; k++ ){
      for ( int j = 0; j < my+1; j++ ){
        for ( int i = 0; i < mx+1; i++ ){
          rootNode[n] = new Node;
          rootNode[n]->key = n;
          if( k != 0 ) {
            if( j != 0 ) {
              if( i != 0  ) rootNode[n]->nbrcell[0] = root[(i-1)+(j-1)*mx+(k-1)*mx*my];
              else          rootNode[n]->nbrcell[0] = NIL;
              if( i != mx ) rootNode[n]->nbrcell[1] = root[(i+0)+(j-1)*mx+(k-1)*mx*my];
              else          rootNode[n]->nbrcell[1] = NIL;
            }
            if( j != my ) {
              if( i != 0  ) rootNode[n]->nbrcell[2] = root[(i-1)+(j+0)*mx+(k-1)*mx*my];
              else          rootNode[n]->nbrcell[2] = NIL;
              if( i != mx ) rootNode[n]->nbrcell[3] = root[(i+0)+(j+0)*mx+(k-1)*mx*my];
              else          rootNode[n]->nbrcell[3] = NIL;
            }
          }
          if( k != mz ) {
            if( j != 0 ) {
              if( i != 0  ) rootNode[n]->nbrcell[4] = root[(i-1)+(j-1)*mx+(k+0)*mx*my];
              else          rootNode[n]->nbrcell[4] = NIL;
              if( i != mx ) rootNode[n]->nbrcell[5] = root[(i+0)+(j-1)*mx+(k+0)*mx*my];
              else          rootNode[n]->nbrcell[5] = NIL;
            }
            if( j != my ) {
              if( i != 0  ) rootNode[n]->nbrcell[6] = root[(i-1)+(j+0)*mx+(k+0)*mx*my];
              else          rootNode[n]->nbrcell[6] = NIL;
              if( i != mx ) rootNode[n]->nbrcell[7] = root[(i+0)+(j+0)*mx+(k+0)*mx*my];
              else          rootNode[n]->nbrcell[7] = NIL;
            }
          }
          n++;
        }
      }
    }

    /** Neighboring Node from cell */
    for ( int in = 0; in < (mx+1)*(my+1)*(mz+1); in++ ){
      for ( int ie = 0; ie < 8; ie ++ ){
        if( rootNode[in]->nbrcell[ie] != NIL ) {
          rootNode[in]->nbrcell[ie]->nbrnode[7-ie] = rootNode[in];
          rootNode[in]->xg[0] = rootNode[in]->nbrcell[ie]->xg[0][7-ie];
          rootNode[in]->xg[1] = rootNode[in]->nbrcell[ie]->xg[1][7-ie];
          rootNode[in]->xg[2] = rootNode[in]->nbrcell[ie]->xg[2][7-ie];
        }
      }
    }

    return key;

  }
  
  std::tuple<Element*,int> zordering( Element* node, Element* nodefront, int noe ){

    if( !node->ibnd ){
      noe++;
      node->key = noe;
      amr::cr8_cnx(node);
      if( nodefront == NIL ) {
        front = node;
        node->prev = NIL;
      }
      else{
        node->prev = nodefront;
        node->prev->next = node;
      }
      nodefront = node;
    }

    return std::forward_as_tuple(nodefront,noe);
  }

  std::tuple<Element*,int> zordering_bnd( Element* node, Element* nodefront, int noe ){

    if( node->ibnd ){
      noe++;
      node->key = noe;
      amr::cr8_cnx(node);
      node->prev = nodefront;
      node->prev->next = node;
      nodefront = node;
    }

    return std::forward_as_tuple(nodefront,noe);
  }
  
}
