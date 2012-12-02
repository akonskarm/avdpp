//    Copright (C) 2012, Alexandros Konstantinakis-Karmis
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA,
//    or download the License terms from prep.ai.mit.edu/pub/gnu/COPYING-2.0.
//
//    Contact:
//    --------
//    Alexandros Konstantinakis-Karmis
//    ΕρΓΑ Laboratory of Algebraic and Geometric Algorithms
//    Division of Theoretical Computer Science
//    Department of Informatics and Telecommunications 
//    National and Kapodistrian University of Athens
//    akonskarm@gmail.com

#include "wspd.h"

//Find out the number of concentric circles that must be drawn for a particular wspd_center/distance
//combination. AMM proposes a particular formula, but we must make sure that our concentric circles
//are always fully contained inside the unit hypercube
//This can be done by following this analysis. max(center) denotes the maximum coordinate of center
//2^i * l + max(center) >= 1 => i >= log2((1-max(center))/l) => i = ceil... floor the last valid
//min(center) - 2^i * l <= 0 => i >= log2(min(center)/l) => i = ceil... floor the last valid
//TODO: make a struct and compute this while i fill the midpoint
int calc_circs_num(ANNpoint center, ANNcoord distance) {
  ANNcoord min_elem = 1.0;
  ANNcoord max_elem = 0.0;
  for (int i=0; i<Globals::dim; i++) {
    if (min_elem > center[i])
      min_elem = center[i];
    if (max_elem < center[i])
      max_elem = center[i];
  }
  
  ANNcoord lower = floor(log2(min_elem/distance));
  ANNcoord upper = floor(log2((1-max_elem)/distance));
  ANNcoord amm   = ceil(log2(Globals::c1/Globals::epsilon));

  return std::min(lower,std::min(upper,amm));
}

void find_larger_smaller(FSTnode *a, FSTnode *b) {
  if (a->radius > b->radius) {
    test_separation(a,b);
  } else {
    test_separation(b,a);
  }
}

void insert_new_pair(FSTnode *larger, FSTnode *smaller, ANNcoord distance) {
  ANNpoint midpoint = annAllocPt(Globals::dim);
  for (int i=0; i<Globals::dim; i++) {
    midpoint[i] = (larger->center[i] + smaller->center[i])/2.0;
  }
  ANNpoint vector = annAllocPt(Globals::dim);
  for (int i=0; i<Globals::dim; i++) {
    vector[i] = larger->center[i] - smaller->center[i];
  }
  Globals::wspd_centers.push_back(midpoint);
  Globals::wspd_vectors.push_back(vector);
  Globals::wspd_distances.push_back(distance);
  Globals::wspd_circs_num.push_back(calc_circs_num(midpoint,distance));
}

//so I want to examine the two non-empty children of a node
//2 leaves are always well separated
//but leaf radius = 0!
void test_separation(FSTnode *larger, FSTnode *smaller) {
  assert ((larger->radius > smaller->radius) || (larger->radius == 0 && smaller->radius == 0));
  ANNcoord distance = sqrt(annDist(Globals::dim, larger->center, smaller->center));
  ANNcoord separation_factor = 4 * larger->radius;
  if (distance > separation_factor) {
    insert_new_pair(larger,smaller,distance);
  } else {
    find_larger_smaller(larger->left, smaller);
    find_larger_smaller(larger->right, smaller);
  }
}

void find_nodes(FSTnode *node) {
  if (node->radius!=0) {
    //INNER NODE
    find_larger_smaller(node->left, node->right);
    find_nodes(node->left);
    find_nodes(node->right);
  }
}

void find_wsp(FSTnode *root) {
  find_nodes(root);
}
