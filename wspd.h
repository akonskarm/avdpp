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

#ifndef ANN
#define ANN
#include <ANN/ANN.h>
#endif
#ifndef GLOBALS
#define GLOBALS
#include "globals.h"
#endif
#ifndef FST
#define FST
#include "fst.h"
#endif

//Find out the number of concentric circles that must be drawn for a particular wspd_center/distance
//combination. AMM proposes a particular formula, but we must make sure that our concentric circles
//are always fully contained inside the unit hypercube
//This can be done by following this analysis. max(center) denotes the maximum coordinate of center
//2^i * l + max(center) >= 1 => i >= log2((1-max(center))/l) => i = ceil... floor the last valid
//min(center) - 2^i * l <= 0 => i >= log2(min(center)/l) => i = ceil... floor the last valid
//TODO: make a struct and compute this while i fill the midpoint
int calc_circs_num(ANNpoint center, ANNcoord distance);

void find_larger_smaller(FSTnode *a, FSTnode *b);

void insert_new_pair(FSTnode *larger, FSTnode *smaller, ANNcoord distance);

//so I want to examine the two non-empty children of a node
//2 leaves are always well separated
//but leaf radius = 0!
void test_separation(FSTnode *larger, FSTnode *smaller);

void find_nodes(FSTnode *node);

void find_wsp(FSTnode *root);
