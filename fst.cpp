//Fair Split Tree
//Help structure to compute the WSPD
//Based on the slow, sequential algorithm described in the Callahan-Kosaraju paper
#include "fst.h"
  FSTnode::FSTnode(FSTnode *l, FSTnode *r, std::vector<ANNpoint> points): left(l),right(r)
  {  
    
    //minimum enclosing circle code
    Miniball mb(Globals::dim);
    for (unsigned int i=0; i<points.size(); i++) {
      mb.check_in(points[i]);
    }
    mb.build();
    center = mb.center();// annAllocPt(Globals::dim);
    radius = sqrt(mb.squared_radius());
  }
  
  FSTnode::~FSTnode()
  {
    delete left;
    delete right;
    delete[] center;
  }
  
  int FSTnode::is_leaf() {
    return ((left==NULL) && (right==NULL));
  }

FSTnode *build_fst(std::vector<ANNpoint> points) {
  if (points.size() == 1) {
    return new FSTnode(0,0,points);
  } else {
    //for (int i=0; i<points.size(); i++) {
    //  printPt(std::cout, points[i]);
    //}
    //std::cout << "---" <<std::endl;
    //INNER NODE
    //This is simple: take the points, find the max dimension,
    //cut the max dimension in half, divide the points in left and
    //right child. recurse.
    
    //Bounding box is the box that contains all points.
    //ll is the lower left corner and rr is the upper right corner
    ANNpoint bbox_ll = annAllocPt(Globals::dim);
    ANNpoint bbox_rr = annAllocPt(Globals::dim);

    //find the bounding box for all input points by computing the min
    //and max element on every axis
    for (int axis=0; axis<Globals::dim; axis++)
    {
      ANNcoord min_elem = 1.0;
      ANNcoord max_elem = 0.0;
      for (unsigned int i=0; i<points.size(); i++) {
        if (points[i][axis] < min_elem) {
          min_elem = points[i][axis];
        }
        if (points[i][axis] > max_elem) {
          max_elem = points[i][axis];
        }
      }
      bbox_ll[axis] = min_elem;
      bbox_rr[axis] = max_elem;
    }
    
    /*std::cout << " ";
    printPt(std::cout, bbox_ll);
    std::cout << " ";
    printPt(std::cout, bbox_rr);*/
    
    //Find the longest axis
    ANNcoord longest_axis_length = 0.0;
    int longest_axis = 0;
    for (int axis=0; axis<Globals::dim; axis++) {
      if (bbox_rr[axis]-bbox_ll[axis]>longest_axis_length) {
        longest_axis_length = bbox_rr[axis]-bbox_ll[axis];
        longest_axis = axis;
      }
    }
   
    //coord of where the bbox is cut
    ANNcoord half_axis = (bbox_ll[longest_axis]+bbox_rr[longest_axis])/2.0;
    //std::cout << longest_axis << "   " << half_axis << std::endl;
    //creating two new point sets for left and right child...READABILITY
    std::vector<ANNpoint> left_points, right_points;
    for (unsigned int i=0; i<points.size(); i++) {
      if (points[i][longest_axis] <= half_axis) {
        left_points.push_back(points[i]);
      } else {
        right_points.push_back(points[i]);
      }
    }
    
    //for (int i=0; i<left_points.size(); i++) {
    //  printPt(std::cout, left_points[i]);
    //}
    //std::cout << std::endl;
    //for (int i=0; i<right_points.size(); i++) {
    //  printPt(std::cout, right_points[i]);
    //}
    //std::cout << " --- " << std::endl;
    //assert(left_points.size() > 0);
    //assert(right_points.size() > 0);
    //assert(points.size() == left_points.size() + right_points.size());
    delete[] bbox_ll;
    delete[] bbox_rr;
    return new FSTnode(build_fst(left_points), build_fst(right_points),points);
  }
}
