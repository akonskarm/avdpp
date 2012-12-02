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

#include "query.h"

ANNpoint query_pt;
Vc vcin;
clock_t time_kdt_exa=0;
clock_t time_kdt_app=0;
clock_t time_avd=0;
clock_t tmp;
ANNidxArray  nnIdx_exa;
ANNdistArray dists_exa;
ANNidxArray  nnIdx_app;
ANNdistArray dists_app;

double maximum_eps=100;
double epsilon_kdt;
double epsilon_avd;
myset* avd_repr;
ANNcoord ndist;
ANNcoord tmp_dist;
ANNpoint qbcenter;
QNode* sqb;
int siz;
int *error_kdt;
int *error_avd;

myset* findLeaf() {
//void findLeaf(ANNpoint* center) {
  double qb_side=0.25;
  int i;
  for (i=0; i<Globals::dim; i++) {
    //hypercube? :( TODO
    qbcenter[i]=0.5;
  }
  sqb = Globals::qb;

  while (1) {
    if (sqb->children == 0) {
      return sqb->representatives;
//      return;
    } else {
      int child_addr=0;
      int it=1;
      for (i=0; i<Globals::dim; i++) {
        if (query_pt[i] > qbcenter[i]) {
          child_addr = child_addr | it;
          qbcenter[i] += qb_side;
        } else {
          qbcenter[i] -= qb_side;
        }
        it = it << 1;
      }
      sqb=&(sqb->children[0][child_addr]);
      qb_side = qb_side/2.0;
    }
  }
}

void insert_error_kdt(ANNcoord e) {
  if (e < pow(10,-6)) {
    error_kdt[0]+=1;
  } else if (e <= Globals::real_e) {
    error_kdt[int(e*10)+1]+=1;
  } else {
    error_kdt[siz-1]+=1;
  }
}
void insert_error_avd(ANNcoord e) {
  if (e < pow(10,-6)) {
    error_avd[0]+=1;
  } else if (e <= Globals::real_e) {
    error_avd[int(e*10)+1]+=1;
  } else {
    error_avd[siz-1]+=1;
  }
}
void test_point() {
    tmp=clock();
    Globals::kdtree->annkSearch(query_pt, 1, nnIdx_exa, dists_exa, 0.0);
    time_kdt_exa+=clock()-tmp;

    tmp=clock();
    Globals::kdtree->annkSearch(query_pt, 1, nnIdx_app, dists_app, Globals::real_e);
    time_kdt_app+=clock()-tmp;
    
    tmp=clock();
    avd_repr = findLeaf();
    ndist = 10.0;//arbitrary long
    for (myset::iterator it = (*avd_repr).begin(); it != (*avd_repr).end(); it++) {
      tmp_dist = annDist(Globals::dim,query_pt,Globals::points[*it]);
      if (tmp_dist < ndist) {
        ndist=tmp_dist;
      }
    }
    time_avd+=clock()-tmp;
 
    epsilon_kdt = sqrt(dists_app[0]/dists_exa[0])-1;
    epsilon_avd = sqrt(ndist/dists_exa[0])-1;
 
    if (maximum_eps > epsilon_avd)
      maximum_eps = epsilon_avd;

    insert_error_kdt(epsilon_kdt);
    insert_error_avd(epsilon_avd);
}

void asdf(QNode* qb, ANNpoint* qb_center, ANNcoord qb_side) {
  if (qb->children == 0) {
    //trying to get about 10^9 queries in total
    for (int j=0; j<pow(10,5)/Globals::ov_leaves; j++) {
      for (int i=0; i<Globals::dim; i++) {
        //TODO: random in the area defined by the box
        query_pt[i] = (*qb_center)[i] + rand()/double(RAND_MAX)*qb_side - qb_side/2.0;
        //TODO: totally random
        //query_pt[i] = rand()/double(RAND_MAX);
    }
    test_point();
    }
  } else {
    int it;
    ANNpointArray child_centers = annAllocPts(Globals::pow2dim, Globals::dim);
    for (int i=0; i<Globals::pow2dim; i++) {
      it = 1;
      for (int j=0; j<Globals::dim; j++) {
        if ((i&it) != 0) {
          child_centers[i][j] = (*qb_center)[j]+(qb_side/4.0);
        } else {
          child_centers[i][j] = (*qb_center)[j]-(qb_side/4.0);
        }   
        it = it << 1;
      }   
      //new_children.push_back(find_all_reprs(qb->children[0] + i,&child_centers[i],side/2.0));
      asdf(&(qb->children[0][i]),child_centers+i,qb_side/2.0);
    }
    annDeallocPts(child_centers);
  }
}

void cart_product(Vc& in) {
    int idx;
    Vd vd;
    // Start all of the iterators at the beginning.
    for (int i=0; i<Globals::dim; i++) {
      Digits d = {in.begin(), in.end(), in.begin()};
      vd.push_back(d);
      query_pt[i]=*(in.begin());
    }

    while(1) {
        // Construct your first product vector by pulling 
        // out the element of each vector via the iterator.
        //Vc result;
        //for(Vd::const_iterator it = vd.begin();
        //    it != vd.end();
        //    it++) {
        //    result.push_back(*(it->me));
        //}
        //out.push_back(result);
        //printPt(std::cout,query_pt);
        test_point();
        // Increment the rightmost one, and repeat.

        // When you reach the end, reset that one to the beginning and
        // increment the next-to-last one. You can get the "next-to-last"
        // iterator by pulling it out of the neighboring element in your
        // vector of iterators.
        idx=0;
        for(Vd::iterator it = vd.begin(); ; ) {
            // okay, I started at the left instead. sue me
            ++(it->me);
            if(it->me == it->end) {
                if(it+1 == vd.end()) {
                    // I'm the last digit, and I'm about to roll
                    return;
                } else {
                    // cascade
                    it->me = it->begin;
                    query_pt[idx]=*(it->me);
                    ++it;
                }
            } else {
                // normal
                query_pt[idx]=*(it->me);
                break;
            }
            ++idx;
        }
    }
}

void query() {

  if (Globals::side_multiplier == 1) {
    Globals::real_e = 0.5;
  } else {
    Globals::real_e = 1.0;
  }
//  Globals::real_e = Globals::side_multiplier * Globals::epsilon;

//  std::cout << "Real_e: " << Globals::real_e << std::endl;

//  siz=(1+log2(Globals::side_multiplier))*5+2;
  siz=int(Globals::real_e*10) + 2;
  error_kdt = new int[siz];
  error_avd = new int[siz];
  for (int i=0; i<siz; i++) {
    error_kdt[i]=0;
    error_avd[i]=0;
  }

  std::cout << std::endl;

  qbcenter = annAllocPt(Globals::dim);
  query_pt = annAllocPt(Globals::dim);

  nnIdx_exa = new ANNidx[1];
  dists_exa = new ANNdist[1];
  nnIdx_app = new ANNidx[1];
  dists_app = new ANNdist[1];

/*  double stuff = 0.1;
  for (int i=0; i<pow(1000000,1.0/Globals::dim); i++) {
    vcin.push_back(stuff);
    stuff+=0.8/pow(1000000,1.0/Globals::dim);
  }
 
  cart_product(vcin);*/

//  for (double num = 0.0; num < 1.0; num+=0.1) {
//    vcin.push_back(num);
//  }
//  cart_product(vcin);

  srand(time(0));

//  for (int i=0; i<pow(10,1); i++) {
//    for (int j=0; j<Globals::dim; j++) {
//      queryPt[j] = (rand()/double(RAND_MAX));
//    }

//  }
 
/*  std::cout << std::endl; 
  std::cout << "TIME EXA: " << time_kdt_exa << std::endl;
  std::cout << "TIME APP: " << time_kdt_app << std::endl;
  std::cout << "TIME AVD: " << time_avd << std::endl;*/
  asdf(Globals::qb, &Globals::hypercube_center, 1.0);


//  std::cout << "\tKD tree Exact" << "\tKD tree Approx" << "\tAVD" << std::endl;
  std::cout << time_kdt_exa << " " << time_kdt_app << " " << time_avd << " "; //<< " " << CLOCKS_PER_SEC << " ";

//  std::cout << "\tKDT\tAVD" << std::endl;
  for (int i=0; i<siz; i++) {
    std::cout << error_kdt[i] << " ";
    error_kdt[i]=0;
  }
  for (int i=0; i<siz; i++) {
    std::cout << error_avd[i] << " ";
    error_avd[i]=0;
  }
//  time_kdt_exa=0;
//  time_kdt_app=0;
//  time_avd=0;


//  vcin.clear();
//  stuff = 0.4;
//  for (int i=0; i<pow(1000000,1.0/Globals::dim); i++) {
//    vcin.push_back(stuff);
//    stuff+=0.2/pow(1000000,1.0/Globals::dim);
//  }

//  cart_product(vcin);
//  std::cout << time_kdt_exa << " " << time_kdt_app << " " << time_avd << " ";

//  for (int i=0; i<7; i++) {
//    std::cout << error_kdt[i] << " ";
//    error_kdt[i]=0;
//  }
//  for (int i=0; i<7; i++) {
//    std::cout << error_avd[i] << " ";
//    error_avd[i]=0;
//  }


  delete[] nnIdx_exa;
  delete[] dists_exa;
  delete[] nnIdx_app;
  delete[] dists_app;
  delete[] error_kdt;
  delete[] error_avd;

  annDeallocPt(qbcenter);
  annDeallocPt(query_pt);
}
