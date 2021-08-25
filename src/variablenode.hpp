#ifndef VARIABLENODE_HPP
#define VARIABLENODE_HPP

class Element;

class Node {
  public:
    int key;
    int noc;
    Node *prev, *next;
    /** The location*/
    double xg[3];
    /** The pointer of a neiboring cell */
    Element *nbrcell[8];
    /** velocity */
    double uo[3];
    /** Pressure */
    double p=0e0, pp=0e0; 
    /** Property */
    double den; 
    double gamm, visc, vist;
};
extern Node *NILNode, **rootNode;
extern Node *frontnode;

#endif
