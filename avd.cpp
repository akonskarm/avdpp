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

#include <iostream>
#include <cmath>
//#include <vector>
#include <fstream>

#include "avd.h"
#include "scale.h"
#include "fst.h"
#include "wspd.h"
#include "overlap.h"
#include "query.h"
#include <sys/resource.h>

void check_mem() {
  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret;
  ret = getrusage(who, &usage);
  if (ret == 0) {
    std::cout << "   " << usage.ru_maxrss << std::endl;
    std::cout << "   " << usage.ru_ixrss << std::endl;
    std::cout << "   " << usage.ru_idrss << std::endl;
    std::cout << "   " << usage.ru_isrss << std::endl;
  }
}

void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

void getArgs(int argc, char **argv) {
  if (argc <= 1) {
    std::cerr << "Usage:\n\n" << "  avd [-d dimensions] [-n points] [-r representatives] [-t threads] [-s side multiplier]";
    std::cerr << "      \n\n" << "  -d: number of dimensions, at least 3";
    std::cerr << "      \n\n" << "  -n: number of points, at least 2";
    std::cerr << "      \n\n" << "  -r: number of representatives per leaf, at least 1";
    std::cerr << "      \n\n" << "  -t: number of threads, at least 1";
    std::cerr << "      \n\n" << "  -s: side_multiplier, 1 for vanilla squares\n";
    std::cerr << "      \n\n" << "  -f: fill leaves with reprs completely\n";
    std::cerr << "      \n\n" << "  -q: just construct, no queries\n";
    exit(0);
  }
  int i = 1;
  while (i < argc) {
    if (!strcmp(argv[i], "-d")) {
      Globals::dim = atoi(argv[++i]);
      Globals::dimm2 = Globals::dim-2;
      Globals::pow2dim = int(pow(2,Globals::dim));
      Globals::pow2dimm2 = int(pow(2,Globals::dimm2));
      Globals::pow2dimm2 = int(pow(2,Globals::dimm2));
      Globals::c2 = 20.0 * Globals::dim;
    } else if (!strcmp(argv[i], "-n")) {
      Globals::len = atoi(argv[++i]);
    } else if (!strcmp(argv[i], "-r")) {
      Globals::reprs_num = atoi(argv[++i]);
    } else if (!strcmp(argv[i], "-t")) {
      Globals::threads_num = atoi(argv[++i]);
    } else if (!strcmp(argv[i], "-s")) {
      Globals::side_multiplier = atoi(argv[++i]);
    } else if (!strcmp(argv[i], "-q")) {
      Globals::queries = false;
    } else if (!strcmp(argv[i], "-f")) {
      Globals::fill = true;
    } else {
      std::cerr << "Unrecognized option.\n";
      exit(1);
    }
    i++;
  }
  //hack lol
  if ((unsigned int) Globals::len <= Globals::reprs_num) {
    std::cout << "Number of points must be greater than number of representatives" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if ((Globals::dim < 3) || (Globals::len < 2) || (Globals::reprs_num < 1) || (Globals::threads_num < 1) || (Globals::side_multiplier < 1)) {
    std::cout << "Wrong parameter" << std::endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char **argv) {
//  std::cout << "SIZE " << sizeof(boost::mutex) << std::endl;
  //setup();
  double vm, rss;
  getArgs(argc, argv);
  scale();
  //for (int i =0; i<Globals::len; i++) {
  //  printPt(std::cout, Globals::points[i]);
  //}
  
  FSTnode *root = build_fst(Globals::points);
  find_wsp(root);
//  find_wsp(build_fst(Globals::points));

  Globals::hypercube_center = annAllocPt(Globals::dim);
  for (int i=0; i<Globals::dim; i++) {
    Globals::hypercube_center[i]=0.5;
  }

  //make the kdtree
  Globals::ann_points = annAllocPts(Globals::len, Globals::dim);
  
  
  for (int i=0; i<Globals::len; i++) {
    for (int j=0; j<Globals::dim; j++) {
      Globals::ann_points[i][j] = Globals::points[i][j];
    }
  }
  Globals::kdtree = new ANNkd_tree(Globals::ann_points, Globals::len, Globals::dim);

  double timer = GetTimeMs();
  overlap();
  
//  std::cout << "fill is " << Globals::fill << std::endl;

  double stop_timer = GetTimeMs()-timer;
  
  process_mem_usage(vm, rss);
  Globals::ov_nodes=1;
  Globals::ov_leaves=0;
  rec_count(Globals::qb);
  check_tree(Globals::qb,0);

//  std::cout << "Dim" << "\tPoints" << "\tReprs" << "\tSide." << "\tWSPairs" << "\tNodesQT" << "\tLeavesQT" << "\tTotalTime" << "\tMergingTime" << "\tMaxDepth" << std::endl;

  std::cout << Globals::dim << " " << Globals::len << " " << Globals::reprs_num << " " << 
  Globals::side_multiplier << " " << Globals::wspd_centers.size() << " " << Globals::ov_nodes << " " << 
//  Globals::ov_leaves << " " << Globals::final_set.size() << " " << (int) stop_timer/1000 << " " << (int) Globals::merging_time/1000 << " " << (int) rss/1000 << " " << 
  Globals::ov_leaves << " " << (int) stop_timer/1000 << " " << (int) Globals::merging_time/1000 << " " << 
  //Globals::sets[0].size() << " " << Globals::max_depth << " ";
  Globals::max_depth << " ";

  if (Globals::queries)
    query();

  annDeallocPts(Globals::ann_points);
  delete Globals::kdtree;
  annClose();

  delete[] Globals::nnIdx;
  delete[] Globals::dists;

  delete root;
  delete[] Globals::move;
  delete[] Globals::hypercube_center;
  //TODO
  for (unsigned int i=0; i<Globals::wspd_centers.size(); i++) {
    delete[] Globals::wspd_centers[i];
  }
  for (int i=0; i<Globals::len; i++) {
    annDeallocPt(Globals::points[i]);
  }

/*  std::map<myset,myset*>::iterator it;
  for (int i=0; i<Globals::threads_num; i++) {
    for (it=Globals::sets[i].begin(); it!=Globals::sets[i].end(); it++) {
      delete (*it).second;
    }
  }
  delete[] Globals::sets;*/
  
  delete Globals::qb;
  delete[] Globals::qbs;

  std::cout << " " << std::endl;

  return EXIT_SUCCESS;
}
