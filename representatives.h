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

/* representatives.h
 * Purpose here is to to a depth first in the quadtree and fill all leaves
 * with a number of representatives. Secondly, if a node's children are all
 * leaves, we try to merge them together and reduce that path's length by one.
 * We use both QLeaf_tmp and QLeaf. In both of them children is always null.
 * Their difference is that representatives is an std::set in QLeaf_tmp
 * and a pointer to an array of ints for QLeaf. std::set is convenient because
 * every element is unique, but a pointer to an array is the most succint.
 */

/* Both needed by ANN. Only ever going to ask for 1 nearest neighbor.
 * Remember that we save the indexes for representatives etc
 */


#include <ANN/ANN.h>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>


#include "globals.h"
#ifndef REPRESENTATIVES_H
#define REPRESENTATIVES_H

//ANNidxArray nnIdx = new ANNidx[1]; //index of nearest neighbor
//ANNdistArray dists = new ANNdist[1]; //square distance of query and NN
/* Used when we try to merge leaves  */
extern __thread int leaves_count;
//extern __thread myset union_reprs;
extern int merges;
extern int nothings;
extern int normals;
extern int leaves;
extern int inner;
/* Key and value contain the same information, but the int** is more succint
 * This map allows us to keep every set only once and let leaves point to them
 */

void depth_all_reprs(QNode* qb);
myset* find_all_reprs(QNode* qb, ANNpoint* center, double side,  ANNkd_tree* kdtree, ANNmin_k* asdf,int id);
myset* fill_in_reprs(QNode* qb, ANNpoint* center, double side,  ANNkd_tree* kdtree, ANNmin_k* asdf);
void replace_reprs(QNode* qb, ANNpoint* center, double side,ANNidxArray nnIdx, ANNdistArray dists);

void fill_vectors(QNode* qb, ANNpoint* center, double side);

void representatives();

struct Slave {
  void operator()();
  Slave(int id);
  int id;
  ANNmin_k* asdf;
  ANNkd_tree* kdtree;
  myset* fil_in_reprs(QNode* qb, ANNpoint* center, double side);
};

void dothereprs();

#endif
