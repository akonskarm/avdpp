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

/* Overlap, the meat of the process
      Globals::pow2dim = int(pow(2,Globals::dim));
 * A generalized Bressenham algorithm that aims to find the overlapping squares
 * The idea for 2D is that, given a number of concentric circles, I'm going to 
 * draw the circumference of a circle and fill it and then draw the circumference
 * of 2 circles per ring and fill their difference, where "fill" refers to the act
 * of finding the center of each square that overlaps the circle or ring.
 * These squares are inserted into a quadtree, our final search structure.
 *
 * For dimensions bigger than 2D, we generalize thusly:
 * Given that the equation for a hypersphere is 
 * x^2 + y^2 + (a_0^2 + a_1^2 + ... + a_d^2) = r^2
 * we just need to do what we did in 2D for each combination of the possible values
 * of a_0^2,...,a_d^2 for a radius of sqrt(r^2 - (a_0^2+...+a_d^2)
 * We can also notice that negative and positive values have the same square and take
 * advantage of this mirroring to do less processing
 * The algorithm is as follows:
 * 1) Generate all the possible numbers for a particular circle or ring. That is either
      0...radius with step side, radius to 2radius with step 2side etc. This is done in the MAIN LOOP
   2) Produce the cartesian product of 1) for (Globals::dim - 2) instances
   3) For each combination, compute the new radius (sqrt(r^2 - sum)). Some combinations are valid, others are not
      Steps 2 and 3 are in cart_product
   4) Find the new "center" coordinates by adding and subtracting the value from the center coordinates.
      So if the center is (5,5) and I get (0,1) as relative coords, it makes the center (5,6).
      Do the cartesian product of all these tuples and take the unique
   5) For every such new center, do a bressenham and fill it
   6) for every bressenham square, insert it!
 */

#include <iterator>
#include <iostream>
#include <algorithm>
#include <sys/types.h>
#include "representatives.h"
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/tss.hpp>
#include <ANN/ANN.h>
#include "globals.h"
#include <vector>

#include "ann_1.1.2/src/pr_queue_k.h"

#ifndef VC
#define VC
typedef std::vector<ANNcoord> Vc;
struct Digits {
    Vc::const_iterator begin;
    Vc::const_iterator end;
    Vc::const_iterator me;
};
typedef std::vector<Digits> Vd;
#endif

#ifndef OVERLAP_H
#define OVERLAP_H
typedef std::vector<Vc> Vvc;
 extern int couter;
extern boost::mutex io_mutex;
//extern boost::mutex leaf_mutex;
//extern boost::mutex repr_mutex;
//extern boost::thread_specific_ptr<int> ptr;
extern ANNcoord denominator;
extern boost::mutex gi_mutex;
extern int global_index;
 
//Discretization function for the side of each square of the grid
ANNcoord side_discr(ANNcoord num, ANNcoord extra);

//Discretization function for the coordinates of each circle or ring
ANNcoord center_discr(ANNcoord num, ANNcoord mult, ANNcoord extra);
 
//Discretization function for the radius
ANNcoord radius_discr(ANNcoord num, ANNcoord mult);

//Helper functions, nothing to see here
void printVc(Vc v);

void printVvc(Vvc vv);


//ANNdist
ANNcoord sum_sqr(int dim, Vc& p);

struct Worker
{
  int cart_prod_counter;
  //std::map<myset,myset*> sets;
  bool core;
  Vc all_values, radiuses, radiuses_sqr, sides, pos_neg_result, result;
  Vvc pos_neg, pos_neg_final;
  int pos_neg_final_size;
  ANNpoint qb_center;
  long long ov_leaves,ov_nodes;
  Vd vd1,vd2;
  std::vector<ANNpoint> centers;
  ANNpointArray leaf_centers;
  int* leaf_reprs;
  ANNkd_tree* kdtree;
  ANNidxArray nnIdx;
  ANNdistArray dists;
  QNode* qb;
  ANNmin_k* asdf;

  int id;
  int private_index;
  
  //returns a unique next index or -1
  //should be another class with a single object
  //(singleton?) with private global_index... oh well!
  int get_index();
  
  Worker(int id);

  void find_leaf2(const ANNpoint* p, const ANNcoord side, const int repr );

  void cart_product(const Vc& in, Vc& result, const int circ);
  //http://stackoverflow.com/questions/5279051/how-can-i-create-cartesian-product-of-vector-of-vectors
  void raster_circle_wrapper(const Vc &result, const int rc, const ANNcoord sum_of_squares);
  /*Standard bressenham, hacked to fill the circle and draw rings etc*/
  void raster_circle(const ANNcoord x0, const ANNcoord y0, const ANNcoord side, const ANNcoord sum_squares, const int circ);
  void find_leaf_wrapper(const ANNcoord x, const ANNcoord y, const ANNcoord side);
  //int representative;
  void find_leaf(const ANNcoord side,const int idx);
  void operator()();
};
 
void overlap();
#endif
