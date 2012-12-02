//Fair Split Tree
//Help structure to compute the WSPD
//Based on the slow, sequential algorithm described in the Callahan-Kosaraju paper


#ifndef FST_H
#define FST_H
#include <ANN/ANN.h>
#include "globals.h"
#include <vector>
#include "Miniball_dynamic_d.h"



class FSTnode
{
public:
  FSTnode *left;
  FSTnode *right;
  //Minimum enclosing Circle of points of the subtree defined
  //with this node as its root
  ANNpoint center;
  ANNcoord radius;
  
  FSTnode(FSTnode *l, FSTnode *r, std::vector<ANNpoint> points);
  
  ~FSTnode();
  
  int is_leaf();
};

FSTnode *build_fst(std::vector<ANNpoint> points);
#endif
