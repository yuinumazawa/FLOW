#ifndef VARIABLE_HPP
#define VARIABLE_HPP
#include <fstream> 

#define ELEMENT

extern int dim;
extern int dimn;
extern int dimb; // 1->001, 3->011, 7->111 
extern int mx, my, mz;
extern double xtot, ytot, ztot;
extern double delt;
extern double reso;
extern std::FILE *file;
extern double *crs_val, *x, *b;
extern int n, *crs_col, *crs_index, noc;
extern double reso;

class Node;

class Element {
  public:
    /** key */
    int key;
    /** The indenter of cell*/
    int xi[3];
    /** The size of cell*/
    double xl[3];
    /** The location of cell center*/
    double xc[3];
    /** The location of cell edges*/
    double xg[3][8];
    /** The volume of cell*/
    double vol;
    /** The surface area of the neiboring cell*/
    double sa[6];
    /** The pointer of parent and children */
    Element *parent, *child[8];
    int childid;
    /** The pointer of a neiboring cell with the same level */
    Element *nbr[6];
    /** The pointer of a neiboring cell*/
    Element *cnx[6];
    /** The pointer of a neiboring node */
    Node *nbrnode[8];
    /** The distance to the neiboring cell*/
    double dis_cnx[6];
    /** The number of the neiboring cells*/
    int num_cnx[6];
    /** My level on the hierarchy of AMR*/
    int lvl;
    Element *prev, *next;
    /** refine number */
    int  rfn = 0;
    /** Flag of object*/
    bool obj = false;
    /** Flag of leaf*/
    bool leaf = true;
    /** value*/
    double val, valo;
    /** coeficient for CFD*/
    double ap, apo, sp, sm;
    double an[6];
    /** boundary */
    bool ibnd = false;
    bool wall = false;
    bool ibndn[6];
    bool outside[6];
    bool movingwalln[6];
    bool walln[6];
    bool inflown[6];
    bool outflown[6];
    bool slipn[6];
    /** velocity */
    double u[3], uo[3]; // amr
    double uf[3][6];
    double dudx[3][3];
    double dufdx[3][3][6];
    double divu;
    double rotu;
    /** Pressure */
    double p=0e0, pp=0e0; //amr
    double pf[6], ppf[6];
    double dPdx[3], dPPdx[3];
    double res = 0e0;
    /** Property */
    double den, deno; // kg/m3 amr
    double gamm, visc, vist; // amr
};

extern Element *NIL, **root, **ele, *front;

#endif
