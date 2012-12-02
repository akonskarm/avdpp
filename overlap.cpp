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

#include "overlap.h"
boost::mutex io_mutex;
ANNcoord denominator = Globals::c2 * Globals::gamma; 
boost::mutex gi_mutex;
int global_index=-1;
int couter = 0;
//Discretization function for the side of each square of the grid
ANNcoord side_discr(ANNcoord num, ANNcoord extra) {
  return pow(2, ceil(log2(num) + extra));
}

//Discretization function for the coordinates of each circle or ring
ANNcoord center_discr(ANNcoord num, ANNcoord mult, ANNcoord extra) {
  return mult * (floor(num / mult) + extra);
}
 
//Discretization function for the radius
ANNcoord radius_discr(ANNcoord num, ANNcoord mult) {
  return mult * round(num / mult);
}

//Helper functions, nothing to see here
void printVc(Vc v) {
  std::cout << "(";
  
  for (unsigned int i=0; i<v.size()-1; i++) {
    std::cout << v[i] << ", ";
  }
  std::cout << v[v.size()-1] << ")" << std::endl;
}

void printVvc(Vvc vv) {
  for (unsigned int i=0; i<vv.size(); i++) {
    printVc(vv[i]);
  }
  std::cout << std::endl;
  std::cout << std::endl;
}

//ANNdist
ANNcoord sum_sqr(int dim, Vc& p)
{
  register int d;
  register ANNcoord dist;

  dist=0;
  for (d=0; d<dim; d++)  {
    dist+=p[d]*p[d];//pow(p[d],2);
  }
  return dist;
}

  
  //returns a unique next index or -1
  //should be another class with a single object
  //(singleton?) with private global_index... oh well!
  int Worker::get_index() {
	boost::mutex::scoped_lock lock(gi_mutex);
	global_index++;
	if (global_index < (int) Globals::wspd_centers.size()) {
	  return global_index;
	} else {
	  return -1;
	}
  }
  
  Worker::Worker(int id) : id(id) {
    //TODO:delete
    asdf = new ANNmin_k(1);
    //nnIdx = new ANNidx[1];
    //dists = new ANNdist[1];
    leaf_centers=annAllocPts(pow(2,Globals::dimm2),Globals::dim);
    leaf_reprs=new int[Globals::pow2dimm2];
    qb_center = annAllocPt(Globals::dim);
    std::vector<int> indices;
    //leaf_center=annAllocPt(Globals::dim);
    for (int i=0; i<Globals::dimm2; i++) {
      indices.push_back(0);
      result.push_back(0.0);
      pos_neg_result.push_back(0.0);
      Vc pos_neg_dim;
      pos_neg_dim.push_back(0.0);
      pos_neg_dim.push_back(0.0);
      pos_neg.push_back(pos_neg_dim);
    }
    for (int i=0; i<pow(2,Globals::dimm2); i++) {
      Vc pos_neg_dim_final;
      for (int j=0; j<Globals::dim; j++)
        pos_neg_dim_final.push_back(0.0);
      pos_neg_final.push_back(pos_neg_dim_final);
    }
    for (int i=0; i<Globals::dimm2; i++) {
      Digits d;
      vd2.push_back(d);
    }
    for (int i=0; i<Globals::dimm2; i++) {
      Digits d;
      vd1.push_back(d);
    }

    kdtree = new ANNkd_tree(Globals::ann_points, Globals::len, Globals::dim);
  }


void Worker::cart_product(const Vc& in, Vc& result, const int circ) {
  /* I have an array of numbers from 0 to some radius saved in "in"
   * At first I need the cartesian product of these values
   * for Globals::dimm2 = Globals::dim - 2 instances
  */

/*  for (unsigned int i=0; i<result.size(); i++) {
    result[i]=0.0;
  }*/
  int idx;
  //the indices
  //squares_sum
  ANNcoord ss;

  //Start all of the iterators at the beginning
  for (int i=0; i<Globals::dimm2; i++) {
    vd1[i].begin=in.begin();
    vd1[i].end=in.end();
    vd1[i].me=in.begin();
  }
  while(1) {
    // Increment the rightmost one, and repeat.
    // When you reach the end, reset that one to the beginning and
    // increment the next-to-last one. You can get the "next-to-last"
    // iterator by pulling it out of the neighboring element in your
    // vector of iterators.
    idx = 0;
    for(Vd::iterator it = vd1.begin(); ; ) {
      // okay, I started at the left instead. sue me
      ++(it->me);
      if (it->me == it->end) {
        //not needed all the time obviously
        if(it + 1 == vd1.end()) {
          // I'm the last digit, and I'm about to roll
          return;
        } else {
          // cascade
          it->me = it->begin;
          result[idx] = *(it->me);
          ++it;
        }
      } else {
        // normal
        result[idx] = *(it->me);
        break;
      }
      ++idx;
    }

    //if the current cartesian product is valid for a new sphere of radius circ:
    ss=sum_sqr(Globals::dimm2, result);
    if (ss <= radiuses_sqr[circ]) {
      raster_circle_wrapper(result, circ, ss);
      //depth_all_reprs(Globals::qbs[id]);
    }

    
    cart_prod_counter++;

    if (cart_prod_counter >= (int) in.size()) {
      ov_nodes=1;
      ov_leaves=0;
      rec_count(Globals::qbs[id]);
//      std::cout << "in here"<< std::endl;
      if (ov_nodes+ov_leaves > pow(2,20)) {
        //find_all_reprs(Globals::qbs[id], &Globals::hypercube_center, 1.0, kdtree, asdf, id);
        depth_all_reprs(Globals::qbs[id]);
    //    std::cout << "find_all_reprs" << std::endl;
      }
      cart_prod_counter=0;
    }
  }
}

//http://stackoverflow.com/questions/5279051/how-can-i-create-cartesian-product-of-vector-of-vectors
void Worker::raster_circle_wrapper(const Vc &result, const int rc, const ANNcoord sum_of_squares) {
  //2.do all the possible combinations of the extra dims...
  for (int j=0; j<Globals::dimm2; j++) {
    pos_neg[j][0] = centers[rc][j+2] + result[j];
    pos_neg[j][1] = centers[rc][j+2] - result[j];
  }

  //3. for each combination of these extra dims...
  int idx=0;
  int id_comb=0;
  int id_dimm=0;
  for (Vvc::const_iterator it = pos_neg.begin(); it != pos_neg.end(); ++it) {
    vd2[idx].begin=(*it).begin();
    vd2[idx].end=(*it).end();
    vd2[idx].me=(*it).begin();
    pos_neg_result[idx]=*(it->begin());
    idx++;
  }
  while(1) {
    for (id_dimm=0; id_dimm<Globals::dimm2; id_dimm++)
      //pos_neg_final[id_comb][id_dimm+2] = pos_neg_result[id_dimm];
      leaf_centers[id_comb][id_dimm+2] = pos_neg_result[id_dimm];
    id_comb++;

    idx=0;
    for(Vd::iterator it = vd2.begin(); ; ) {
      ++(it->me);
      if (it->me == it->end) {
        if(it+1 == vd2.end()) {
          
          //sort(pos_neg_final.begin(), pos_neg_final.end(), vcsort);
          //vvcend = std::unique(pos_neg_final.begin(),pos_neg_final.end(),vccomp);
          //pos_neg_final_size = vvcend - pos_neg_final.begin();
          //printVvc(pos_neg_final);
          //pos_neg_final
          //pos_neg_final_size = pos_neg_final.size();
          raster_circle(centers[rc][0],centers[rc][1],sides[rc],sum_of_squares,rc);
          return;
        } else {
          it->me = it->begin;
          pos_neg_result[idx] = *(it->me);
          ++it;
        }
      } else {
        pos_neg_result[idx] = *(it->me);
        break;
      }
    ++idx;
    }
  }
}

/*Standard bressenham, hacked to fill the circle and draw rings etc*/
void Worker::raster_circle(const ANNcoord x0, const ANNcoord y0, const ANNcoord side, const ANNcoord sum_squares, const int circ) {
   ANNcoord radius = sqrt(radiuses_sqr[circ]-sum_squares);
   ANNcoord radius_prev = 0;
   ANNcoord border;
   if (core) {
     find_leaf_wrapper(x0 ,y0 , side);
     border=0.0;
   } else {
     radius_prev = sqrt(radiuses_sqr[circ-1] - sum_squares);
     border= ceil(radius_prev/side)*side;
   }

   ANNcoord f = side - radius;
   ANNcoord f_prev = side - radius_prev;
   ANNcoord ddF_x = side;
   ANNcoord ddF_y = -2.0*radius;
   ANNcoord ddF_y_prev = -2.0*radius_prev; //or simply -radius... d'oh
   ANNcoord x = 0.0;
   ANNcoord y = radius;

   for (ANNcoord i = side + border; i <= radius; i += side) {
     find_leaf_wrapper(x0, y0 + i, side);
     find_leaf_wrapper(x0, y0 - i, side);
     find_leaf_wrapper(x0 + i, y0, side);
     find_leaf_wrapper(x0 - i, y0, side);
   }
   
   while (true) {
     if (f >= 0) {
       y -= side;
       ddF_y += 2.0 * side;
       f += ddF_y;
     }
     if (!core) {
       if (x<radius/2.0) {
         if (f_prev >= 0) {
           border -= side;
           ddF_y_prev += 2.0 * side;
           f_prev += ddF_y_prev;
         }
       } else {
         border=0.0;
       }
     }
     x += side;
     ddF_x += 2.0 * side;
     f += ddF_x;
     if (!core) {
       f_prev += ddF_x;
     }
     if(x>=y) {
       break;
     }

     for (ANNcoord i = side + border; i < y+side; i += side) {
       find_leaf_wrapper(x0 + x, y0 + i, side);
       find_leaf_wrapper(x0 + x, y0 - i, side);
       find_leaf_wrapper(x0 - x, y0 + i, side);
       find_leaf_wrapper(x0 - x ,y0 - i, side);
     }
     for (ANNcoord i=y; i>sqrt(2)/2*radius; i-=side) {
       find_leaf_wrapper(x0 + i, y0 + x, side);
       find_leaf_wrapper(x0 + i, y0 - x, side);
       find_leaf_wrapper(x0 - i, y0 + x, side);
       find_leaf_wrapper(x0 - i, y0 - x, side);
     }

   }
   if (fabs(x-y) < pow(1,-5) * x) {
     find_leaf_wrapper(x0 + x, y0 + y, side);
     find_leaf_wrapper(x0 + x, y0 - y, side);
     find_leaf_wrapper(x0 - x, y0 + y, side);
     find_leaf_wrapper(x0 - x, y0 - y, side);
   }
}

void Worker::find_leaf_wrapper(const ANNcoord x, const ANNcoord y, const ANNcoord side) {
  //Just a little dumb function to fill in the missing information of x and y
  //Call the inserter to quadtree for each "true" set of coordinates"

  for (int i=0; i<Globals::pow2dimm2; i++) {
    leaf_centers[i][0]=x;
    leaf_centers[i][1]=y;
  }
 
 
//  for (int i=0; i<Globals::pow2dimm2; i++) {
//    kdtree->annkSearch(leaf_centers[i],1,nnIdx,dists);
//    leaf_reprs[i]=nnIdx[0];

//    leaf_reprs[i]=kdtree->ann1Search(leaf_centers[i],asdf);
//  }


//  {boost::mutex::scoped_lock lock(leaf_mutex);
  for (int i=0; i<Globals::pow2dimm2; i++) {
    if (leaf_centers[i][0]*Globals::wspd_vectors[private_index][0] + 
        leaf_centers[i][1]*Globals::wspd_vectors[private_index][1] + 
     //   leaf_centers[i][3]*Globals::wspd_vectors[private_index][3] + 
        leaf_centers[i][2]*Globals::wspd_vectors[private_index][2] < 2*side)
    find_leaf(side,i);
//  }}
  }
}

//int representative;
//or replace idx with a &leaf_centers thing TODO
void Worker::find_leaf(const ANNcoord side, const int idx) {
  couter++;
  register int i;
//  int representative = -1;
  qb=Globals::qbs[id];
  for (i=0; i<Globals::dim; i++) {
    //or hypercube center? mmm TODO
    qb_center[i]=0.5;
  }
  register double qb_side = 0.5;//to megethos twn paidiwn tou current node

  //why scoped_lock
  register int child_addr;
  register int it;
  while (qb_side >= side) {
//    boost::mutex::scoped_lock lock(qb->node_mutex);
    //qb->node_mutex.lock();
    if (qb->children == 0) {
      if (qb->depth == 0) {
//        std::cout << "DEPTH ZERO" << std::endl;
        qb->bear_children();
      } else {
//        std::cout << "DEPTH NOT ZERO" << std::endl;
//        std::cout << (int) log2(2*qb_side/side)  << " " << qb->depth << std::endl;
        //open the subtree, insert the square, merge
        if ((int) log2(2*qb_side/side) > qb->depth) {
//          std::cout << "LOG" << std::endl;
          
          qb->bear_children();
          for (int i=0; i<Globals::pow2dim; i++) {
            qb->children[0][i].depth = qb->depth-1;
          }
          qb->depth=0;
        } else {
          return;
        }
      }
    }
    child_addr=0;
    it=1;
    for (i=0; i<Globals::dim; i++) {
      if (leaf_centers[idx][i] > qb_center[i]) {
        child_addr = child_addr | it;
        qb_center[i] += qb_side;
      } else {
        qb_center[i] -= qb_side;
      }
      it = it << 1;
    }
    qb=&(qb->children[0][child_addr]);
    qb_side = qb_side/2.0;
  }

}

  void Worker::operator()() {
    for (private_index=get_index(); private_index!=-1; private_index=get_index()) {
      {boost::mutex::scoped_lock lock(io_mutex);
      std::cout << "THREAD " << id << " GOT " << private_index+1 << "/" << Globals::wspd_centers.size() << std::endl;}
    //Find the side length, radius, radius squared and center of circle/ring for
    //each one of the circle/rings that I'm going to draw
      for (int j=0; j<=Globals::wspd_circs_num[private_index]; j++) {
        //1 as argument because I insert parents of leaves
        //This means that when I do find a square with this side, I'm going to call
        //bear_children method again. This essentially reduces processing time by 3/4
        sides.push_back(side_discr(pow(2,j)*Globals::wspd_distances[private_index]/denominator,3));//1
//        sides.push_back(side_discr(pow(2,j)*Globals::wspd_distances[private_index]/denominator,0));
        //make radius an exact multiple of side
        radiuses.push_back(radius_discr(pow(2,j)*Globals::wspd_distances[private_index], sides[j]));
        radiuses_sqr.push_back(radiuses[j]*radiuses[j]);
        //Center of circle should fall neatly on the center of a grid square
        ANNpoint center = annAllocPt(Globals::dim);
        for (int k=0; k<Globals::dim; k++) {
          center[k] = center_discr(Globals::wspd_centers[private_index][k], sides[j], 0.5);
        }
        centers.push_back(center);
      }

      //first iteration is for the core, the inner circle which is completely filled
      //the rest of the iteration are rings, essentially 2 concentric cirles where I fill
      //the difference
      core=true;
      for (int circ=0; circ<=Globals::wspd_circs_num[private_index]; circ++) {
        //Generate all the possible values of a coordinate, based on side and radius
        //for example if radius is 8 and side is 2 this is going to be 0,2,4,6,8
        for (ANNcoord v=0.0; v<=radiuses[circ]; v+=sides[circ]) {
          all_values.push_back(v);
        }
        //cart_prod_counter=0;
        cart_product(all_values, result, circ);
        if (Globals::qbs[id]->children!=0)
        depth_all_reprs(Globals::qbs[id]);

//        find_all_reprs(Globals::qbs[id], &Globals::hypercube_center, 1.0, kdtree, asdf, id);
        
        //CLEANUP
        all_values.clear();
        core=false;
      }
    
      //CLEANUP
      radiuses.clear();
      radiuses_sqr.clear();
      sides.clear();
      for (unsigned int d=0; d<centers.size(); d++)
        annDeallocPt(centers[d]);
      centers.clear();

    }
    std::cout << "THREAD " << id << " DIES" << std::endl;
    annDeallocPt(qb_center);
    annDeallocPts(leaf_centers);
    delete[] leaf_reprs;
    pos_neg_result.clear();
    pos_neg_final.clear();
    vd1.clear();
    vd2.clear();
    delete asdf;
    delete kdtree;
//    delete qb;
    return;
  }

//TODO: lock to map???
//TODO: ann1Search sto reprs
/*void merge_nodes(QNode *ag, QNode *it) {
  if ((ag->children != 0) && (it->children !=0)) {
    for (int i=0; i<Globals::pow2dim; i++) {
      merge_nodes(ag->children[0] + i, it->children[0] + i);
    }
  } else if ((ag->children == 0) && (it->children !=0)) {
    if (ag->representatives != 0) {
      ag->remove_representatives();
    }
    ag->bear_children();
    for (int i=0; i<Globals::pow2dim; i++) {
      merge_nodes(ag->children[0] + i, it->children[0] + i);
    }
  } else if ((ag->children==0) && (it->children==0)) {
    myset union_reprs;
    if (ag->representatives!=0) {
      for (myset::iterator i = ag->representatives[0].begin(); i != ag->representatives[0].end(); i++) {
        union_reprs.insert(*i);
      }
    }
    if (it->representatives!=0) {
      for (myset::iterator i = it->representatives[0].begin(); i != it->representatives[0].end(); i++) {
        union_reprs.insert(*i);
      }
    }

    if (Globals::sets[0].count(union_reprs) == 0) {
      ag->representatives = new myset(union_reprs);
      Globals::sets[0][union_reprs] = ag->representatives;
    } else {
      ag->representatives = Globals::sets[0][union_reprs];
    }
    //std::cout << union_reprs.size() << " " << ag->representatives[0].size() << std::endl;
  }
}*/
void merge_nodes(QNode *ag, QNode *it) {
  if ((ag->depth != 0) && (it->depth != 0)) {
    ag->depth = std::max(ag->depth,it->depth);
  } else if ((ag->depth == 0) && (it->depth != 0)) {
    if (ag->children != 0) {
      std::cout << "error merging?" <<std::endl;
    }
    ag->depth = it->depth; //MERGE with previous
  } else if ((ag->children==0) && (it->children==0)) {
    for (int i=0; i<Globals::pow2dim; i++) {
      merge_nodes(ag->children[0]+i,it->children[0]+i);
    }
  }
}
void make_it_deeper(QNode *qb) {
  if (qb->children == 0) {
    qb->depth+=4;
  } else {
    for (int i=0; i<Globals::pow2dim; i++) {
      make_it_deeper(qb->children[0]+i);
    }
  }
}

/*void replace_reprs(QNode* qb, ANNidxArray* nnIdx, ANNdistArray* dists) {
  if (qb->children == 0) {
    if (qb->representatives[0].size() < Globals::reprs_num) {
      Globals::kdtree->annkSearch(query_pt, reprs_num, nnIdx_exa, dists_exa, 0.0);
      
    }

    if (Globals::final_set.count(qb->representatives[0]) == 0) {
      Globals::final_set[qb->representatives[0]] = qb->representatives;
    } else {
      qb->representatives = Globals::final_set[qb->representatives[0]];
    }
  } else {
    for (int i=0; i<Globals::pow2dim; i++) {
      replace_reprs(qb->children[0]+i);
    }
  }
}*/


void overlap() { 
  Globals::qbs = new QNode*[Globals::threads_num];
  Globals::sets = new std::map<myset,myset*>[Globals::threads_num];
  for (int i = 0; i<Globals::threads_num; i++) {
    Globals::qbs[i] = new QNode();
  }

  denominator = Globals::c2 * Globals::gamma; 
  denominator /= Globals::side_multiplier;
  
  boost::thread* thr[Globals::threads_num];
  for (int i = 0; i < Globals::threads_num; ++i)
    thr[i] = new boost::thread(Worker(i));
  for (int i = 0; i < Globals::threads_num; ++i) {
    thr[i]->join();
    delete thr[i];
  }
 
//  for (int i=0; i<Globals::threads_num; i++) {
//    rec_count(Globals::qbs[i], 1.0);
//    std::cout << Globals::ov_nodes << " " << Globals::ov_leaves << std::endl;
//    Globals::ov_nodes=1;
//    Globals::ov_leaves=0;
//  }
 
  double merging_timer = GetTimeMs();
  for (int i=1; i<Globals::threads_num; i++) {
    std::cout << "Merging " << i << std::endl;
    merge_nodes(Globals::qbs[0],Globals::qbs[i]);
    delete Globals::qbs[i];
  }
  Globals::qb = Globals::qbs[0];
  //depth_all_reprs(Globals::qb);
  make_it_deeper(Globals::qb);

  //TODO
  //std::cout << "start filling" << std::endl;

//  std::cout << "before fill vectors" << std::endl;
//  fill_vectors(Globals::qb, &Globals::hypercube_center, 1.0);

/*  for (int i=0; i<10; i++) {
    printPt(std::cout,*Globals::qt_centers[i]);
    std::cout<<std::endl;
    std::cout << Globals::qt_sides[i] << std::endl;
  }*/

  Globals::ov_nodes=1;
  Globals::ov_leaves=0;
  rec_count(Globals::qb);
  //TODO
  ////std::cout << " ovnodes " << Globals::ov_nodes << " " << Globals::ov_leaves << std::endl;
  ////std::cout << "LEAF: " << couter << std::endl;

  dothereprs();

ANNmin_k* asdf = new ANNmin_k(1);

//  fill_in_reprs(Globals::qb, &Globals::hypercube_center, 1.0, Globals::kdtree, asdf);
  find_all_reprs(Globals::qbs[0], &Globals::hypercube_center, 1.0, Globals::kdtree, asdf,0);

  Globals::merging_time = GetTimeMs() - merging_timer;

  ANNidxArray  nnIdx = new ANNidx[Globals::reprs_num];
  ANNdistArray dists = new ANNdist[Globals::reprs_num];
  //replace_reprs(Globals::qb, &Globals::hypercube_center, 1.0, nnIdx, dists); CHANGE
  delete[] nnIdx;
  delete[] dists; 
//TODO
//  std::cout << "done replacing" << std::endl; 
  /*std::map<myset,myset*>::iterator it;
  for (int i=0; i<Globals::threads_num; i++) {
    for (it=Globals::sets[i].begin(); it!=Globals::sets[i].end(); it++) {
      delete (*it).second;
    }
  }*/
  delete[] Globals::sets;

  //representatives();
}
