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
#include "globals.h"
#include <cmath>
#ifndef QUERY_H
#define QUERY_

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


extern ANNpoint query_pt;
extern Vc vcin;
extern clock_t time_kdt_exa;
extern clock_t time_kdt_app;
extern clock_t time_avd;
extern clock_t tmp;
extern ANNidxArray  nnIdx_exa;
extern ANNdistArray dists_exa;
extern ANNidxArray  nnIdx_app;
extern ANNdistArray dists_app;

extern double epsilon_kdt;
extern double epsilon_avd;
extern myset* avd_repr;
extern ANNcoord ndist;
extern ANNcoord tmp_dist;
extern ANNpoint qbcenter;
extern QNode* sqb;
extern int siz;
extern int* error_kdt;
extern int* error_avd;

myset* findLeaf();

void insert_error_kdt(ANNcoord e);
void insert_error_avd(ANNcoord e);
void test_point();
void cart_product(Vc& in);


void asdf(QNode* qb, ANNpoint* qb_center, ANNcoord qb_side);
void query();
#endif
