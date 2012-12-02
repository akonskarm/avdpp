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

#include "representatives.h"

#include <algorithm>
__thread int leaves_count;

boost::mutex thread_mutex;
boost::mutex acc_mutex;
int thread_index=-1;

/* Used when we try to merge leaves  */
/* Key and value contain the same information, but the int** is more succint
 * This map allows us to keep every set only once and let leaves point to them
 */

//this is the merger for leaves with no reprs
void depth_all_reprs(QNode* qb) {
  bool all_children_leaves=true;
  //id==depth
  bool all_children_same_id=true;
  for (int i=0; i<Globals::pow2dim; i++) {
    if (qb->children[0][i].children != 0) {
      all_children_leaves=false;
      break;
    } else {
      if (qb->children[0][0].depth != qb->children[0][i].depth) {
        all_children_same_id = false;
      }
    }
  }

  if (all_children_leaves) {
    if (all_children_same_id) {//same depth!
      qb->depth = qb->children[0][0].depth + 1;
      qb->remove_children();
    } else {
      //inner node forevah
      qb->depth=0;
    }
  } else {
    for (int i=0; i<Globals::pow2dim; i++) {
      if (qb->children[0][i].children != 0) {
        depth_all_reprs(qb->children[0]+i);
      }
    }
    all_children_leaves=true;
    all_children_same_id=true;
    //copy - paste
    for (int i=0; i<Globals::pow2dim; i++) {
      if (qb->children[0][i].children != 0) {
        all_children_leaves=false;
        break;
      } else {
        if (qb->children[0][0].depth != qb->children[0][i].depth) {
          all_children_same_id = false;
        }
      }
    }
    if ((all_children_leaves) && (all_children_same_id)) {
      qb->depth = qb->children[0][0].depth + 1;
      qb->remove_children();
    } else {
      qb->depth=0;
    }
  }
}

//process all that shit in vectors.... very gay
void fill_vectors(QNode* qb, ANNpoint* center, double side) {
  if (qb->children == 0) {
    Globals::qt_nodes.push_back(qb);
    Globals::qt_centers.push_back(center);
    Globals::qt_sides.push_back(side);
    //std::cout << side << std::endl;
  } else {
    int it;
    ANNpointArray child_centers = annAllocPts(Globals::pow2dim, Globals::dim);
    for (int i=0; i<Globals::pow2dim; i++) {
      it = 1;
      for (int j=0; j<Globals::dim; j++) {
        if ((i&it) != 0) {
          child_centers[i][j] = (*center)[j]+(side/4.0);
        } else {
          child_centers[i][j] = (*center)[j]-(side/4.0);
        }
        it = it << 1;
      }
      fill_vectors(qb->children[0]+i,&child_centers[i],side/2.0);
    }
    annDeallocPts(child_centers);
  }
}

//fill in the reprs in the given QNode
myset* fill_in_reprs(QNode* qb, ANNpoint* center, double side, ANNkd_tree* kdtree, ANNmin_k* asdf) {
  if ((qb->depth == 0) && (qb->children == 0)) {
    //time to fill in the reprs my dear!
    //if (qb->representatives != 0) {
    //  std::cout << "something wicked this way comes" << std::endl;
    //}
    myset rs;
    rs.insert(kdtree->ann1Search(*center,asdf));
    {
    //if (Globals::final_set.count(rs) == 0) {
    if (Globals::final_set.count(rs) == 0) {
      qb->representatives = new myset(rs);
      Globals::final_set[rs] = qb->representatives;
    } else {
      qb->representatives = Globals::final_set[rs];
    }
    }
    return qb->representatives;
  } else {
    if (qb->depth != 0) {
      qb->bear_children();
      for (int i=0; i<Globals::pow2dim; i++) {
        qb->children[0][i].depth = qb->depth-1;
      }
      qb->depth=0;
    }

    myset union_reprs;
    int it;
    //A local container to consume the result of our recursive calls
    //std::vector<Globals::QLeaf*> new_children;
    myset* new_children[Globals::pow2dim];
    //recurse to children
    ANNpointArray child_centers = annAllocPts(Globals::pow2dim, Globals::dim);
    for (int i=0; i<Globals::pow2dim; i++) {
      it = 1;
      for (int j=0; j<Globals::dim; j++) {
        if ((i&it) != 0) {
          child_centers[i][j] = (*center)[j]+(side/4.0);
        } else {
          child_centers[i][j] = (*center)[j]-(side/4.0);
        }
        it = it << 1;
      }
      new_children[i] = fill_in_reprs(qb->children[0]+i,&child_centers[i],side/2.0, kdtree, asdf);
    }
    annDeallocPts(child_centers);
    leaves_count = 0;
    //check every child and insert union all reprs 
    for (int i=0; i<Globals::pow2dim; i++) {
      if (new_children[i] != 0) {
        leaves_count++;
        //std::set doesn't have a [] operator. Instead we use the iterators
        for (myset::iterator it = (*new_children[i]).begin();
             it != (*new_children[i]).end();
             ++it) {
          union_reprs.insert(*it);
        }
      }
    }
    if ((leaves_count == Globals::pow2dim) && (union_reprs.size() <= Globals::reprs_num)) {
      //Create the leaf and transfer the new set of representatives
      if (Globals::final_set.count(union_reprs) == 0) {
        //implies there are no reprs in this node.
        Globals::final_set[union_reprs] = new myset(union_reprs);
      }
      qb->representatives = Globals::final_set[union_reprs];
      qb->remove_children();
      return qb->representatives;
    } else {
      return 0;
    }
  }

}

//Should be a method inside Worker
myset* find_all_reprs(QNode* qb, ANNpoint* center, double side, ANNkd_tree* kdtree, ANNmin_k* asdf, int id) {
  if (qb->children == 0) {
    /* If the current node has no children, it's a leaf. Therefore we need
     * to calculate a number of representatives and return a QLeaf_tmp
     */
    if (qb->representatives == 0) {
      //std::cout << "UNOREPRS" << std::endl;
      //COPY
      myset rs;
      rs.insert(kdtree->ann1Search(*center,asdf));
      if (Globals::sets[id].count(rs) == 0) {
        qb->representatives = new myset(rs);
        Globals::sets[id][rs] = qb->representatives;
      } else {
        qb->representatives = Globals::sets[id][rs];
      }
      //exit(0);
    }
    return qb->representatives;
  } else {
    /* This is an inner node. First we recurse to each leaf.
     * Notice that we need to compute the coordinates of the centers
     * of each child, since that information is not stored in QNode
     */
    myset union_reprs;
    int it;
    //A local container to consume the result of our recursive calls
    //std::vector<Globals::QLeaf*> new_children;
    myset* new_children[Globals::pow2dim];
    //Just a fancy array of arrays. We have pow2dim children and each one
    //needs dim coordinates to express its center
    ANNpointArray child_centers = annAllocPts(Globals::pow2dim, Globals::dim);

    for (int i=0; i<Globals::pow2dim; i++) {
      it = 1;
      for (int j=0; j<Globals::dim; j++) {
        if ((i&it) != 0) {
          child_centers[i][j] = (*center)[j]+(side/4.0);
        } else {
          child_centers[i][j] = (*center)[j]-(side/4.0);
        }
        it = it << 1;
      }
      new_children[i] = find_all_reprs(qb->children[0]+i,&child_centers[i],side/2.0, kdtree, asdf, id);
    }
    annDeallocPts(child_centers);
    //We now have a vector of size pow2dim, full of QNode and QLeaf_tmp
    //Count the leaves and at the same time take the union of their
    //representatives
    leaves_count = 0;
    //union_reprs.clear();
    //check every child and insert union all reprs 
    for (int i=0; i<Globals::pow2dim; i++) {
      if (new_children[i] != 0) {
        leaves_count++;
        //std::set doesn't have a [] operator. Instead we use the iterators
        for (myset::iterator it = (*new_children[i]).begin();
             it != (*new_children[i]).end();
             ++it) {
          union_reprs.insert(*it);
        }
      }
    }

    //if all children are leaves and the union of all representatives is not
    //bigger than max, I can merge all leaves into one and return that
    if ((leaves_count == Globals::pow2dim) && (union_reprs.size() <= Globals::reprs_num)) {
      //Create the leaf and transfer the new set of representatives
//      if (qb->representatives !=0) {
//      std::cout << "halp" << std::endl;
//      }
      if (Globals::sets[id].count(union_reprs) == 0) {
        //implies there are no reprs in this node.
        Globals::sets[id][union_reprs] = new myset(union_reprs);
      }
      qb->representatives = Globals::sets[id][union_reprs];
      qb->remove_children();
      return qb->representatives;
    } else {
      return 0;
    }
  }
}

//Should be a method inside Worker
//this one both fills out the leaves and also merges all reprs in one map
void replace_reprs(QNode* qb, ANNpoint* center, double side, ANNidxArray nnIdx, ANNdistArray dists) {
  if (qb->children == 0) {
    myset rs(qb->representatives[0]);


    if (Globals::fill) {
      if (rs.size() < Globals::reprs_num) {
        Globals::kdtree->annkSearch(center[0], Globals::reprs_num, nnIdx, dists, 0.0);
        for (unsigned int i=0; i<Globals::reprs_num; i++) {
          if ((!repr_exists(&rs, nnIdx[i])) && (rs.size()<Globals::reprs_num)) {
            rs.insert(nnIdx[i]);
          }
        }
      }
    }

    //if (qb->representatives==0) {
    //  std::cout << "err" << std::endl;
    //}
    //This really doesn't work as expected mmm

/*    qb->bear_children();
    for (int i=0; i<Globals::pow2dim; i++) {
      qb->children[0][i].representatives = qb->representatives;
    }
    qb->remove_representatives();
    int it;
    int new_repr;
    ANNpoint child_center = annAllocPt(Globals::dim);
    for (int i=0; i<Globals::pow2dim; i++) {
      it = 1;
      for (int j=0; j<Globals::dim; j++) {
        if ((i&it) != 0) {
          child_center[j] = (*center)[j]+(side/4.0);
        } else {
          child_center[j] = (*center)[j]-(side/4.0);
        }
        it = it << 1;
      }
      
      //new_repr=kdtree->ann1Search(child_center,asdf);
      //new_repr = Globals::kdtree->annkSearch(child_center, 1, nnIdx, dists, 0.0);
      Globals::kdtree->annkSearch(child_center, 1, nnIdx, dists, 0.0);

      new_repr=nnIdx[0];

      myset rs(qb->children[0][i].representatives[0]);
      if (rs.size() == Globals::reprs_num-1) {
        //remove the furthest away
        double mdist=0.0;
        double tmp_dist;
        myset::iterator it_rem;
        for (myset::iterator it = rs.begin(); it != rs.end(); it++) {
          tmp_dist = annDist(Globals::dim,child_center,Globals::points[*it]);
          if (tmp_dist > mdist) {
            mdist=tmp_dist;
            it_rem=it;
          }
        }
        rs.erase(it_rem);
      }
      rs.insert(new_repr);
      if (Globals::final_set.count(rs) == 0)
        Globals::final_set[rs] = new myset(rs);
      qb->children[0][i].representatives = Globals::final_set[rs];
    }*/
    
    if (Globals::final_set.count(rs) == 0)
      Globals::final_set[rs]=new myset(rs);
    qb->representatives = Globals::final_set[rs];

    /*if (Globals::final_set.count(qb->representatives[0]) == 0)
      Globals::final_set[qb->representatives[0]] = new myset(qb->representatives[0]);
    qb->representatives = Globals::final_set[qb->representatives[0]];*/

   // annDeallocPt(child_center);
  } else {
    /* This is an inner node. First we recurse to each leaf.
     * Notice that we need to compute the coordinates of the centers
     * of each child, since that information is not stored in QNode
     */
    int it;
    //A local container to consume the result of our recursive calls
    //std::vector<Globals::QLeaf*> new_children;
    //Just a fancy array of arrays. We have pow2dim children and each one
    //needs dim coordinates to express its center
    ANNpoint child_center = annAllocPt(Globals::dim);

    for (int i=0; i<Globals::pow2dim; i++) {
      it = 1;
      for (int j=0; j<Globals::dim; j++) {
        if ((i&it) != 0) {
          child_center[j] = (*center)[j]+(side/4.0);
        } else {
          child_center[j] = (*center)[j]-(side/4.0);
        }
        it = it << 1;
      }
      replace_reprs(qb->children[0]+i, &child_center, side/2.0, nnIdx, dists);
    }
    annDeallocPt(child_center);
  }
}

//void representatives() {
//  find_all_reprs(Globals::qb, &Globals::hypercube_center, 1.0, kd);
//}


  Slave::Slave(int id) : id(id) {
    asdf = new ANNmin_k(1);
    kdtree = new ANNkd_tree(Globals::ann_points, Globals::len, Globals::dim);
  }


/*  int Slave::get_index(int id) : id(id) {
    boost::mutex::scoped_lock lock(thread_mutex);
    thread_index++;
    if (thread_index < (int) Globals::qt_centers.size()) {
      return thread_index;
    } else {
      return -1;
    }
  }*/


  myset* Slave::fil_in_reprs(QNode* qb, ANNpoint* center, double side) {
    if ((qb->depth == 0) && (qb->children == 0)) {
      //time to fill in the reprs my dear!
      //if (qb->representatives != 0) {
      //}
      myset rs;
      rs.insert(kdtree->ann1Search(*center,asdf));

//      {boost::mutex::scoped_lock lock(thread_mutex);
      if (Globals::final_set.count(rs) == 0) {
        qb->representatives = new myset(rs);
        Globals::final_set[rs] = qb->representatives;
      } else {
        qb->representatives = Globals::final_set[rs];
      }
//      }
      return qb->representatives;
    } else {
      if (qb->depth != 0) {
        qb->bear_children();
        for (int i=0; i<Globals::pow2dim; i++) {
          qb->children[0][i].depth = qb->depth-1;
        }
        qb->depth=0;
      }
  
      myset union_reprs;
      int it;
      //A local container to consume the result of our recursive calls
      //std::vector<Globals::QLeaf*> new_children;
      myset* new_children[Globals::pow2dim];
      //recurse to children
      ANNpointArray child_centers = annAllocPts(Globals::pow2dim, Globals::dim);
      for (int i=0; i<Globals::pow2dim; i++) {
        it = 1;
        for (int j=0; j<Globals::dim; j++) {
          if ((i&it) != 0) {
            child_centers[i][j] = (*center)[j]+(side/4.0);
          } else {
            child_centers[i][j] = (*center)[j]-(side/4.0);
          }
          it = it << 1;
        }
        new_children[i] = fil_in_reprs(qb->children[0]+i,&child_centers[i],side/2.0);
      }
      annDeallocPts(child_centers);
      leaves_count = 0;
      //check every child and insert union all reprs 
      for (int i=0; i<Globals::pow2dim; i++) {
        if (new_children[i] != 0) {
          leaves_count++;
          //std::set doesn't have a [] operator. Instead we use the iterators
          for (myset::iterator it = (*new_children[i]).begin();
               it != (*new_children[i]).end();
               ++it) {
            union_reprs.insert(*it);
          }
        }
      }
      if ((leaves_count == Globals::pow2dim) && (union_reprs.size() <= Globals::reprs_num)) {
        //Create the leaf and transfer the new set of representatives

        //{boost::mutex::scoped_lock lock(thread_mutex);
        if (Globals::final_set.count(union_reprs) == 0) {
          //implies there are no reprs in this node.
          Globals::final_set[union_reprs] = new myset(union_reprs);
        }
        qb->representatives = Globals::final_set[union_reprs];
        //}
        qb->remove_children();
        return qb->representatives;
      } else {
        return 0;
      }
    }
  }
  
void Slave::operator()() {
  ANNpointArray child_centers = annAllocPts(Globals::pow2dim, Globals::dim);
  int it;
  for (int i=0; i<Globals::pow2dim; i++) {
    it = 1;
    for (int j=0; j<Globals::dim; j++) {
      if ((i&it) != 0) {
        child_centers[i][j] = 0.75;
      } else {
        child_centers[i][j] = 0.25;
      }
      it = it << 1;
    }
    if (i%Globals::threads_num == id) {
      //TODO make this dynamically allocate pieces to threads
      std::cout << "Thread " << id << " gets branch " << i << std::endl;
      fil_in_reprs(Globals::qb->children[0]+i,&child_centers[i],0.5);
    }
  }
  annDeallocPts(child_centers);
}

void dothereprs() {
  boost::thread* thr[Globals::threads_num];
  for (int i = 0; i < Globals::threads_num; ++i)
    thr[i] = new boost::thread(Slave(i));
  for (int i = 0; i < Globals::threads_num; ++i) {
    thr[i]->join();
    delete thr[i];
  }
}
