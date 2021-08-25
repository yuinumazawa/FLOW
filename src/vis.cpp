#include <iomanip> 
#include <iostream> 
#include <vector> 
#include "variable.hpp"

namespace ensg
{

  void output( int type, int noe ){
    int tmp, tmpi = 0; float tmpf;
    if( type == 0 || type == 1 || type == 2 ){
      for( int ie = 0; ie < noe; ie++ ) {
        for( int i = 0; i < 8; i++ ){
          tmpf = ele[ie]->xg[type][i];
          fwrite(&tmpf,sizeof(float),1,file);
        }
      }
    }
    else if( type == 3 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpi = 8*ie + 1;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = 8*ie + 2;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = 8*ie + 4;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = 8*ie + 3;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = 8*ie + 5;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = 8*ie + 6;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = 8*ie + 8;
        fwrite(&(tmpi),sizeof(int),1,file);
        tmpi = 8*ie + 7;
        fwrite(&(tmpi),sizeof(int),1,file);
      }
    }
    else if( type == 4 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->u[0];
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 5 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->u[1];
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 6 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->u[2];
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 7 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->visc;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 8 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->den;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 9 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->deno;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 10 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->p;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 11 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->pp;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 12 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->lvl;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 13 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->rfn;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 14 ){
      for( int ie = 0; ie < noe; ie++ ) {
        if( ele[ie]->obj ) tmpf = 1e0;
        else tmpf = 0e0;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 15 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->divu;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 16 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->rotu;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 17 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = ele[ie]->res;
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    else if( type == 18 ){
      for( int ie = 0; ie < noe; ie++ ) {
        tmpf = 0e0;
        for( int i = 0; i < dimn; i++ ){
          if( ele[ie]->ibndn[i] == true ) {
            if( ele[ie]->inflown[i] == true ) tmpf += 1e0;
            else if( ele[ie]->outflown[i] == true ) tmpf += 2e0;
            else if( ele[ie]->walln[i] == true ) tmpf += 3e0;
            else if( ele[ie]->slipn[i] == true ) tmpf += 5e0;
          }
        }
        fwrite(&tmpf,sizeof(float),1,file);
      }
    }
    
  }
  
  void ensg( int step, int noe ){
    
    int n = 0;
    int inttmp ;
    // case
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
    // case
    std::string str4filecase = "data//amr.case";
    // geometry
    std::string str4filegeo  = "data//mesh.geo";
    // variable
    std::vector<std::string> str4filedat;
    str4filedat.push_back("data//velocity.dat");
    str4filedat.push_back("data//viscosity.dat");
    str4filedat.push_back("data//density.dat");
    str4filedat.push_back("data//densityo.dat");
    str4filedat.push_back("data//pp.dat");
    str4filedat.push_back("data//P.dat");
    str4filedat.push_back("data//lvl.dat");
    str4filedat.push_back("data//slvl.dat");
    str4filedat.push_back("data//obj.dat");
    str4filedat.push_back("data//divu.dat");
    str4filedat.push_back("data//rotu.dat");
    str4filedat.push_back("data//res.dat");
    str4filedat.push_back("data//hoge.dat");
    std::ostringstream sout;
    
    // case file 
    std::ofstream outputfile("data//amr.case");
    outputfile<<"FORMAT"<<std::endl;
    outputfile<<"type: ensight gold"<<std::endl;
    outputfile<<" "<<std::endl;
    outputfile<<"GEOMETRY"<<std::endl;
    outputfile<<"model: mesh.geo****"<<std::endl;
    outputfile<<"VARIABLE"<<std::endl;
    outputfile<<"vector per element: 1 Velocity_[m/s] velocity.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 Viscosity_[Pa_s] viscosity.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 Density_[kg/m3] density.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 Density_Old_[kg/m3] densityo.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 Pressure_increment_[-] pp.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 Pressure_[Pa] P.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 AMR_Level   lvl.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 AMR_Setlvl slvl.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 AMR_obj obj.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 divu_[kg/(m3_s)] divu.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 rotu_[1/s] rotu.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 Relative_residual_[-] res.dat****"<<std::endl;
    outputfile<<"scalar per element: 1 hoge hoge.dat****"<<std::endl;
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
      fwrite(&str ,sizeof(str),1,file);
      fwrite(&str1,sizeof(str),1,file);
      fwrite(&str2,sizeof(str),1,file);
      fwrite(&str3,sizeof(str),1,file);
      fwrite(&str4,sizeof(str),1,file);
      fwrite(&str5,sizeof(str),1,file);
      inttmp=1;
      fwrite(&inttmp,sizeof(int),1,file);
      fwrite(&str6,sizeof(str),1,file);
      fwrite(&str7,sizeof(str),1,file);
      // node
      inttmp=8*(noe);
      fwrite(&inttmp,sizeof(int),1,file);
      for( int id = 0; id < 3; id++ ){
        output( id, noe );
      }
      // element
      fwrite(&str8,sizeof(str),1,file);
      inttmp=noe; fwrite(&inttmp,sizeof(int),1,file);
      output( 3, noe );
    fclose(file);
  
    // variable file 
    int i = 0;
    str4filedat[i] = str4filedat[i] + sout.str();
    file = fopen(str4filedat[i].c_str(),"wb");
    // header
    fwrite(&str9,80,1,file);
    fwrite(&str5,80,1,file);
    inttmp=1;
    fwrite(&inttmp,sizeof(int),1,file);
    fwrite(&str8,80,1,file);
    output( 4, noe );
    output( 5, noe );
    output( 6, noe );
    fclose(file);
  
    for( int i = 1; i < str4filedat.size(); i++ ){
      str4filedat[i] = str4filedat[i] + sout.str();
      file = fopen(str4filedat[i].c_str(),"wb");
      // header
      fwrite(&str9,80,1,file);
      fwrite(&str5,80,1,file);
      inttmp=1;
      fwrite(&inttmp,sizeof(int),1,file);
      fwrite(&str8,80,1,file);
      output( 6+i, noe );
      fclose(file);
    }
  
  }
}
