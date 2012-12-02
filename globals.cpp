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

#include "globals.h"

//typedef std::set<int> myset;

//#include <boost/container/flat_set.hpp>
//typedef boost::container::flat_set<int> myset;
namespace Globals {
  int dimm2;
  int pow2dim;
  int pow2dimm2;

  int dim=0;
  int len=0;
  unsigned int reprs_num=0;
  int threads_num=0;
  int side_multiplier=0;

  bool fill=false;
  bool queries=true;

  int ov_leaves=0;
  int ov_nodes=1;
  int re_leaves=0;
  int re_nodes=0;

  std::map<myset,myset*>* sets;
  std::map<myset,myset*> final_set;
  ANNpointArray ann_points;
  std::vector<ANNpoint> points;
  double epsilon = 0.5;

  double real_e;

  //double c1 = 2.0 + epsilon;
  double c1 = 1.0;
  double gamma = 2;
  ANNpoint move = 0;
  long double scale = -1.0;
  std::vector<ANNpoint> wspd_centers;
  std::vector<ANNpoint> wspd_vectors;
  std::vector<ANNcoord> wspd_distances;
  std::vector<int> wspd_circs_num;

  std::vector<QNode*> qt_nodes;
  std::vector<ANNpoint*> qt_centers;
  std::vector<double> qt_sides;

  ANNpoint hypercube_center = 0;
  double c2;// = 20.0 * dim;

  double merging_time;
  //unsigned int reprs_num = int(ceil(1.0/pow(Globals::gamma * Globals::epsilon,

  int circs_num = int(ceil(log2(Globals::c1/Globals::epsilon)));

  ANNkd_tree* kdtree;
  ANNidxArray nnIdx = new ANNidx[1];
  ANNdistArray dists = new ANNdist[1];
  QNode *qb;// = new QNode();
  QNode **qbs;

  int max_depth=0;
};
QNode::QNode() {
  children = 0;
  representatives = 0;
  depth = 0;
}

void QNode::remove_representatives() {;
  representatives=0;
}

void QNode::remove_children() {
  delete[] children[0];
  delete[] children;
  children=0;
}

QNode::~QNode() {
  if (children != 0) {
    remove_children();
    remove_representatives();
  }
}

void QNode::bear_children() {
  children = new QNode*[1];
  children[0] = new QNode[Globals::pow2dim];
}
  
//namespace Globals
//{
  
//}
void printPt(std::ostream &out, ANNpoint p)			// print point
{
  out.precision(32);
	out << "(" << p[0];
	for (int i = 1; i < Globals::dim; i++)
  {
		out << ", " << p[i];
	}
	out << ")\n";
}

//TODO: check depth!
void check_tree(QNode* qb, int depth) {
  if (qb->children == 0) {
    Globals::max_depth = std::max(Globals::max_depth, depth);
    if (qb->representatives == 0) {
      std::cout << "LEAF NO REPRS" << std::endl;
    } else {
      if (qb->representatives[0].find(-1) != qb->representatives[0].end()) {
        std::cout << "LEAF WITH MINUS ONE" << std::endl;
      }
    }
  } else {
    if (qb->representatives != 0) {
      std::cout << "NODE WITH REPRS" << std::endl;
    }// else {
     // std::cout << "NODE" << std::endl;
   // }
    for (int i=0; i<Globals::pow2dim; i++) {
      check_tree(&(qb->children[0][i]), depth+1);
    }
  }
}

void rec_count(QNode* qb){//, double side) {
    //if (qb->children[0][i].children != 0) {
  if (qb->children != 0) {
    Globals::ov_nodes++;
    for (int i=0; i<Globals::pow2dim; i++)
      //rec_count(&(qb->children[0][i]),side/2.0);
      rec_count(&(qb->children[0][i]));
  } else {
 //     std::cout << "LEAF: " << side << std::endl;
    Globals::ov_leaves++;
//      if (qb->representatives != 0) {
//        std::cout << qb->representatives[0][0] << std::endl;
//        Globals::blah[qb->representatives[0][0]]++;
//      }
  }
}

double GetTimeMs() {
  struct timeval tv;
  gettimeofday(&tv, NULL);

  double ret = tv.tv_usec;
  /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
  ret /= 1000;

  /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
  ret += (tv.tv_sec * 1000);

  return ret;
}

bool repr_exists(myset* representatives, int repr) {
  for (myset::iterator it = (*representatives).begin(); it != (*representatives).end(); it++) {
    if (*it==repr) {
      return true;
    }
  }
  return false;
}
