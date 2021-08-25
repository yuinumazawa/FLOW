#include <iomanip> 
#include <iostream> 
#include <vector> 
#include "variable.hpp"
#include "variablenode.hpp"

namespace ensg_node
{
  int zordering_node( int noe ){
    for( int ie = 0; ie < noe; ie++ ) {
      for( int in = 0; in < 8; in++ ) {
        ele[ie]->nbrnode[in]->key = -1;
      }
    }
    for( int ie = 0; ie < noe; ie++ ) {
      for( int in = 0; in < 8; in++ ) {
        ele[ie]->nbrnode[in]->prev = NILNode;
        ele[ie]->nbrnode[in]->next = NILNode;
        ele[ie]->nbrnode[in]->key  = -1;
      }
    }
    int non = 0, nodekey = 0;
    Node *nodefront = NILNode;
    frontnode = NILNode;
    for( int ie = 0; ie < noe; ie++ ) {
      for( int in = 0; in < 8; in++ ) {
        if( ele[ie]->nbrnode[in] != NILNode && ele[ie]->nbrnode[in]->key == -1 ){
          if( nodefront == NILNode ) {
            nodefront = ele[ie]->nbrnode[in];
            nodefront->prev = NILNode;
            nodefront->key = nodekey;
            frontnode = nodefront;
          }
          else{
            ele[ie]->nbrnode[in]->prev = nodefront;
            ele[ie]->nbrnode[in]->prev->next = ele[ie]->nbrnode[in];
            ele[ie]->nbrnode[in]->key = nodekey;
            nodefront = ele[ie]->nbrnode[in];
          }
          nodekey++;
        }
      }
    }
    non = nodekey;
    return non;
  }

  void ele2node( int non ){
    Node *nodefront = frontnode;
    for( int in = 0; in < non; in++ ) {
      double tmp1 = 0e0, tmp2 = 0e0, tmp3 = 0e0, tmp4=0e0, tmp5=0e0;
      double tmp6 = 0e0;
      double totalvol = 0e0;
      int tmpn = 0;
      nodefront->noc = 0e0;
      for( int ie = 0; ie < 4; ie++ ) {
        nodefront->noc++;
        if( nodefront->nbrcell[ie] != NIL && nodefront->nbrcell[3-ie] != NIL){
          tmpn++;
          totalvol += nodefront->nbrcell[3-ie]->vol;
          tmp1     += nodefront->nbrcell[3-ie]->vol * nodefront->nbrcell[ie]->pp;
          tmp2     += nodefront->nbrcell[3-ie]->vol * nodefront->nbrcell[ie]->p;
          tmp3     += nodefront->nbrcell[3-ie]->vol * nodefront->nbrcell[ie]->uo[0];
          tmp4     += nodefront->nbrcell[3-ie]->vol * nodefront->nbrcell[ie]->uo[1];
          tmp5     += nodefront->nbrcell[3-ie]->vol * nodefront->nbrcell[ie]->uo[2];
          tmp6     += nodefront->nbrcell[3-ie]->vol * nodefront->nbrcell[ie]->visc;
        }
      }
      for( int ie = 0; ie < 4; ie++ ) {
        if( nodefront->nbrcell[ie+4] != NIL && nodefront->nbrcell[3-ie+4] != NIL){
          tmpn++;
          totalvol += nodefront->nbrcell[3-ie+4]->vol;
          tmp1     += nodefront->nbrcell[3-ie+4]->vol * nodefront->nbrcell[ie+4]->pp;
          tmp2     += nodefront->nbrcell[3-ie+4]->vol * nodefront->nbrcell[ie+4]->p;
          tmp3     += nodefront->nbrcell[3-ie+4]->vol * nodefront->nbrcell[ie+4]->uo[0];
          tmp4     += nodefront->nbrcell[3-ie+4]->vol * nodefront->nbrcell[ie+4]->uo[1];
          tmp5     += nodefront->nbrcell[3-ie+4]->vol * nodefront->nbrcell[ie+4]->uo[2];
          tmp6     += nodefront->nbrcell[3-ie+4]->vol * nodefront->nbrcell[ie+4]->visc;
        }
      }
      if( tmpn != 0 ){
        nodefront->pp    = tmp1 / totalvol;
        nodefront->p     = tmp2 / totalvol;
        nodefront->uo[0] = tmp3 / totalvol;
        nodefront->uo[1] = tmp4 / totalvol;
        nodefront->uo[2] = tmp5 / totalvol;
        nodefront->visc  = tmp6 / totalvol;
      }
      else{
        nodefront->pp   = 0e0;
        nodefront->p    = 0e0;
        nodefront->uo[0] = 0e0;
        nodefront->uo[1] = 0e0;
        nodefront->uo[2] = 0e0;
        nodefront->visc = 0e0;
      }

      nodefront = nodefront->next;
    }
  }


  void output( int type, int noe, int non ){
    int tmp, tmpi = 0; float tmpf;
    if( type == 0 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpi = tmpnode->key;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpnode = tmpnode->next;
      }
    }
    else if( type == 1 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpf = tmpnode->xg[0];
        fwrite(&tmpf,sizeof(float),1,file);
        tmpnode = tmpnode->next;
      }
    }
    else if( type == 2 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpf = tmpnode->xg[1];
        fwrite(&tmpf,sizeof(float),1,file);
        tmpnode = tmpnode->next;
      }
    }
    else if( type == 3 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpf = tmpnode->xg[2];
        fwrite(&tmpf,sizeof(float),1,file);
        tmpnode = tmpnode->next;
      }
    }
    else if( type == 4 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpi = ele[ie]->key;
        fwrite(&(tmpi),sizeof(int),1,file);
      }
    }
    else if( type == 5 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpi = ele[ie]->nbrnode[0]->key+1;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = ele[ie]->nbrnode[1]->key+1;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = ele[ie]->nbrnode[3]->key+1;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = ele[ie]->nbrnode[2]->key+1;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = ele[ie]->nbrnode[4]->key+1;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = ele[ie]->nbrnode[5]->key+1;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = ele[ie]->nbrnode[7]->key+1;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = ele[ie]->nbrnode[6]->key+1;
        fwrite(&(tmpi),sizeof(int),1,file);
      }
    }
    else if( type == 6 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpf = tmpnode->uo[0];
        fwrite(&tmpf,sizeof(float),1,file);
        tmpnode = tmpnode->next;
      }
      //for( int ie = 0; ie < noe; ie++ ) {
      //  tmpf = ele[ie]->u[0];
      //  fwrite(&(tmpi),sizeof(int),1,file);
      //}
    }
    else if( type == 7 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpf = tmpnode->uo[1];
        fwrite(&tmpf,sizeof(float),1,file);
        tmpnode = tmpnode->next;
      }
      //for( int ie = 0; ie < noe; ie++ ) {
      //  tmpf = ele[ie]->u[1];
      //  fwrite(&(tmpi),sizeof(int),1,file);
      //}
    }
    else if( type == 8 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpf = tmpnode->uo[2];
        fwrite(&tmpf,sizeof(float),1,file);
        tmpnode = tmpnode->next;
      }
      //for( int ie = 0; ie < noe; ie++ ) {
      //  tmpf = ele[ie]->u[2];
      //  fwrite(&(tmpi),sizeof(int),1,file);
      //}
    }
    else if( type == 9 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpf = tmpnode->visc;
        fwrite(&tmpf,sizeof(float),1,file);
        tmpnode = tmpnode->next;
      }
      //for( int ie = 0; ie < noe; ie++ ) {
      //  tmpf = ele[ie]->visc;
      //  fwrite(&(tmpi),sizeof(int),1,file);
      //}
    }
    else if( type == 10 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpf = tmpnode->p;
        fwrite(&tmpf,sizeof(float),1,file);
        tmpnode = tmpnode->next;
      }
      //for( int ie = 0; ie < noe; ie++ ) {
      //  tmpf = ele[ie]->p;
      //  fwrite(&(tmpi),sizeof(int),1,file);
      //}
    }
    else if( type == 11 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpf = tmpnode->pp;
        fwrite(&tmpf,sizeof(float),1,file);
        tmpnode = tmpnode->next;
      }
      //for( int ie = 0; ie < noe; ie++ ) {
      //  tmpf = ele[ie]->pp;
      //  fwrite(&(tmpi),sizeof(int),1,file);
      //}
    }
    else if( type == 12 ){
      Node *tmpnode = frontnode;
      for( int in = 0; in < non; in++ ) {
        tmpf = tmpnode->noc;
        fwrite(&tmpf,sizeof(float),1,file);
        tmpnode = tmpnode->next;
      }
      //for( int ie = 0; ie < noe; ie++ ) {
      //  tmpf = ele[ie]->pp;
      //  fwrite(&(tmpi),sizeof(int),1,file);
      //}
    }
    
  }
  
  void ensg( int step, int noe, int non ){
    
    int n = 0;
    int inttmp ;
    // geometry
    std::string str4filegeo  = "data//mesh.geo";
    char str[80] = "C Binary";
    char str1[80] = "description line 1";
    char str2[80] = "description line 2";
    char str3[80] = "node id off";
    char str4[80] = "element id off";
    char str5[80] = "part";
    char str6[80] = "description line";
    char str7[80] = "coordinates";
    char str8[80] = "hexa8";
    char str9[80] = "Per_elem scalar values for the EnSight Gold geometry example";
    // variable
    std::vector<std::string> str4filedat;
    str4filedat.push_back("data//velocity.dat");
    str4filedat.push_back("data//viscosity.dat");
    str4filedat.push_back("data//pp.dat");
    str4filedat.push_back("data//P.dat");
    str4filedat.push_back("data//noc.dat");
    std::ostringstream sout;
    // case file 
    std::string str4filecase = "data//amr.case";
    std::ofstream outputfile("data//amr.case");
    outputfile<<"FORMAT"<<std::endl;
    outputfile<<"type: ensight gold"<<std::endl;
    outputfile<<" "<<std::endl;
    outputfile<<"GEOMETRY"<<std::endl;
    outputfile<<"model: mesh.geo****"<<std::endl;
    outputfile<<"VARIABLE"<<std::endl;
    outputfile<<"vector per node: 1 Velocity_[m/s] velocity.dat****"<<std::endl;
    outputfile<<"scalar per node: 1 Viscosity_[Pa_s] viscosity.dat****"<<std::endl;
    outputfile<<"scalar per node: 1 Pressure_increment_[-] pp.dat****"<<std::endl;
    outputfile<<"scalar per node: 1 Pressure_[Pa] P.dat****"<<std::endl;
    outputfile<<"scalar per node: 1 NumberOfCells noc.dat****"<<std::endl;
    outputfile<<" "<<std::endl;
    outputfile<<"TIME"<<std::endl;
    outputfile<<"time set: 1"<<std::endl;
    outputfile<<"number of steps: "<< step+1 << std::endl;
    outputfile<<"filename start number: 0"<<std::endl;
    outputfile<<"filename increment: 1"<<std::endl;
    outputfile<<"time values:"<<std::endl;
    for ( int i = 0; i<=step; i++ ) outputfile<<float(i)<<std::endl;
    outputfile.close();
  
    // geometry file 
    sout << std::setfill('0') << std::setw(4) << step;
    str4filegeo = str4filegeo + sout.str();
    file = fopen(str4filegeo.c_str(),"wb");
      // header
      fwrite(&str ,sizeof(str),1,file); // C Binary
      fwrite(&str1,sizeof(str),1,file); // description line 1"
      fwrite(&str2,sizeof(str),1,file); // description line 2"
      fwrite(&str3,sizeof(str),1,file); // node id off";
      fwrite(&str4,sizeof(str),1,file); // element id off";
      fwrite(&str5,sizeof(str),1,file); //part
      inttmp=1;
      fwrite(&inttmp,sizeof(int),1,file);
      fwrite(&str6,sizeof(str),1,file); //description line 1
      fwrite(&str7,sizeof(str),1,file); //coordinates
      inttmp=non;
      fwrite(&inttmp,sizeof(int),1,file); // nn
      for( int id = 1; id < 4; id++ ){
        output( id, noe, non ); // 0: id_n
      }
      fwrite(&str8,sizeof(str),1,file); //hexa8
      inttmp=noe;
      fwrite(&inttmp,sizeof(int),1,file); // ne
      //output( 4, noe, non ); // id_ne
      output( 5, noe, non );
    fclose(file);
  
    // variable file 
    int i = 0;
    str4filedat[i] = str4filedat[i] + sout.str();
    file = fopen(str4filedat[i].c_str(),"wb");
    // header
    fwrite(&str9,sizeof(str),1,file); //description line 1
    fwrite(&str5,80,1,file); // part
    inttmp=1;
    fwrite(&inttmp,sizeof(int),1,file);
    fwrite(&str8,sizeof(str),1,file); //coordinates
    // vector
    output( 6, noe, non );
    output( 7, noe, non );
    output( 8, noe, non );
    fclose(file);
  
    for( int i = 1; i < str4filedat.size(); i++ ){
      str4filedat[i] = str4filedat[i] + sout.str();
      file = fopen(str4filedat[i].c_str(),"wb");
      // header
      fwrite(&str9,sizeof(str),1,file); //description line 1
      fwrite(&str5,80,1,file); // part
      inttmp=1;
      fwrite(&inttmp,sizeof(int),1,file);
      fwrite(&str8,sizeof(str),1,file); //coordinates
      // scholar
      output( 8+i, noe, non );
      fclose(file);
    }
  
  }
}
