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

#include <ANN/ANN.h>
#include <vector>
#include <boost/interprocess/containers/flat_set.hpp>
//#include <boost/container/flat_set.hpp>
#include <set>
#include <map>
#include <sys/time.h>
#include <ctime>

#ifndef GLOBALS_H
#define GLOBALS_H
typedef boost::interprocess::flat_set<int> myset;
//typedef boost::container::flat_set<int> myset;
//typedef std::set<int> myset;

class QNode {
public:
  QNode** children;
  myset* representatives;
  int depth;

  QNode();
  void remove_representatives();
  void remove_children();
  ~QNode();
  void bear_children();
};

namespace Globals {
  extern int dim,dimm2,pow2dim,pow2dimm2,max_threads,len,threads_num,ov_leaves,ov_nodes,re_leaves,re_nodes,side_multiplier;
 
  extern bool fill;
  extern bool queries;

  extern std::map<myset,myset*>* sets;
  extern std::map<myset,myset*> final_set;
  extern ANNpointArray ann_points;
  extern std::vector<ANNpoint> points;
  extern double epsilon;
  extern double real_e;
  extern double gamma;
  extern ANNpoint move;
  extern long double scale;
  extern std::vector<ANNpoint> wspd_centers;
  extern std::vector<ANNpoint> wspd_vectors;
  extern std::vector<ANNcoord> wspd_distances;
  extern std::vector<int> wspd_circs_num;

  extern std::vector<QNode*> qt_nodes;
  extern std::vector<ANNpoint*> qt_centers;
  extern std::vector<double> qt_sides;


  extern ANNpoint hypercube_center;
  extern double c1;
  extern double c2;

  extern double merging_time;
  //unsigned int reprs_num = int(ceil(1.0/pow(Globals::gamma * Globals::epsilon,
  extern unsigned int reprs_num;

  extern int circs_num;

  extern ANNkd_tree* kdtree;
  extern ANNidxArray nnIdx;
  extern ANNdistArray dists;
  extern QNode *qb;

  extern QNode **qbs;

  extern int max_depth;
};

void printPt(std::ostream &out, ANNpoint p);

//TODO: check depth!
void check_tree(QNode* qb, int depth);

//void rec_count(QNode* qb, double side);
void rec_count(QNode* qb);

double GetTimeMs();
bool repr_exists(myset* representatives, int repr);

#endif
