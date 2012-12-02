//    Copright (C) 1999-2006, Bernd Gaertner
//    $Revision: 1.4 $
//    $Date: 2008/02/11 11:50:05 $
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
//    Bernd Gaertner
//    Institute of Theoretical Computer Science 
//    ETH Zuerich
//    CAB G32.2
//    CH-8092 Zuerich, Switzerland
//    http://www.inf.ethz.ch/personal/gaertner

#include "Miniball_dynamic_d.h"

// Functions
// =========
double mb_sqr (double r) {return r*r;}
 
// Class Declarations
// ==================

// smallest enclosing ball of a set of n points in dimension d

// smallest ball with a set of n <= d+1 points *on the boundary*

// point in dimension d
//class Point;


// Class Definitions
// =================

// Miniball_b
// ----------

  // if you want the copy constructor and assignment operator, please
  // properly implement them yourself. The default copy/assignment
  // semantics fails since a Miniball_b object stores pointers to
  // dynamically allocated memory
  //Miniball_b::Miniball_b (const Miniball_b& mb);
  //Miniball_b::Miniball_b& operator=(const Miniball_b& mb);  

  Miniball_b::Miniball_b(int dim) 
    : d(dim) 
  {
    q0 = new double[d];
    z = new double[d+1];
    f = new double[d+1];
    v = new double*[d+1];
    for (int i=0; i<d+1; ++i) v[i] =  new double[d];
    a = new double*[d+1];
    for (int i=0; i<d+1; ++i) a[i] =  new double[d];   
    c = new double*[d+1];
    for (int i=0; i<d+1; ++i) c[i] =  new double[d];
    sqr_r = new double[d+1];
    reset();
  }

  Miniball_b::~Miniball_b() 
  {
    delete[] sqr_r;
    for (int i=0; i<d+1; ++i) delete[] c[i];
    delete[] c;
    for (int i=0; i<d+1; ++i) delete[] a[i];
    delete[] a;
    for (int i=0; i<d+1; ++i) delete[] v[i];
    delete[] v;
    delete[] f;
    delete[] z;
    delete[] q0;
  }
    

// Miniball
// --------
    
  // private methods
// Point (inline)
// --------------

//class Point {
//private:
  //int d; 
  //double* coord;
   
//public:
  //// default
  //Point(int dim)
    //: d (dim), coord(new double[dim])
  //{}

  //~Point ()
  //{
    //delete[] coord;
  //}
   
  //// copy from Point
  //Point (const Point& p)
    //: d (p.dim()), coord(new double[p.dim()])
  //{
    //for (int i=0; i<d; ++i)
      //coord[i] = p.coord[i];
  //}
   
  //// copy from double*
  //Point (int dim, const double* p)
    //: d (dim), coord(new double[dim])
  //{
    //for (int i=0; i<d; ++i)
      //coord[i] = p[i];
  //}
   
  //// assignment
  //Point& operator = (const Point& p)
  //{
    //assert (d == p.dim());
    //if (this != &p)
      //for (int i=0; i<d; ++i)
	//coord[i] = p.coord[i];
    //return *this;
  //}
   
  //// coordinate access
  //double& operator [] (int i)
  //{
    //return coord[i];
  //}
  //const double& operator [] (int i) const
  //{
    //return coord[i];
  //}
  const ANNcoord *begin(ANNpoint p)
  {
    return &p[0];
  }
  const ANNcoord *end(ANNpoint p)
  {
   return &p[Globals::dim];
  }

  // dimension access
  //int dim() const
  //{
  //  return d;
  //}
//};
   

// Class Implementations
// =====================
    
// Miniball
// --------

Miniball::Miniball(int dim) : d(dim), B(dim) {}

void Miniball::check_in (const ANNpoint& p)
{
  L.push_back(p);
}   

void Miniball::build ()
{
  B.reset();
  support_end = L.begin();
  pivot_mb (L.end());
}
   

void Miniball::mtf_mb (It i)
{
  support_end = L.begin();
  if ((B.size())==d+1) return;
  for (It k=L.begin(); k!=i;) {
    It j=k++;
    if (B.excess(*j) > 0) {
      if (B.push(*j)) {
	mtf_mb (j);
	B.pop();
	move_to_front(j);
      }
    }
  }
}
   
void Miniball::move_to_front (It j)
{
  if (support_end == j)
    support_end++;
  L.splice (L.begin(), L, j);
}

void Miniball::pivot_mb (It i)
{
  It t = ++L.begin();
  mtf_mb (t);
  double max_e, old_sqr_r = -1;
  do {
    It pivot;
    max_e = max_excess (t, i, pivot);
    if (max_e > 0) {
      t = support_end;
      if (t==pivot) ++t;
      old_sqr_r = B.squared_radius();
      B.push (*pivot);
      mtf_mb (support_end);
      B.pop();
      move_to_front (pivot);
    }
  } while ((max_e > 0) && (B.squared_radius() > old_sqr_r));
}
   
double Miniball::max_excess (It t, It i, It& pivot) const
{
  const double *c = B.center(), sqr_r = B.squared_radius();
  double e, max_e = 0;
  for (It k=t; k!=i; ++k) {
    //const double *p = (*k).begin();
    const ANNcoord *p = begin(*k);
    e = -sqr_r;
    for (int j=0; j<d; ++j)
      e += mb_sqr(p[j]-c[j]);
    if (e > max_e) {
      max_e = e;
      pivot = k;
    }
  }
  return max_e;
}
   
ANNpoint Miniball::center () const
{
  ANNpoint temp = annAllocPt(Globals::dim);
  for (int i=0; i<Globals::dim; i++) {
    temp[i]=B.center()[i];
  }
  return temp;//Point(d, B.center());
}
   

double Miniball::squared_radius () const
{
  return B.squared_radius();
}
   
   

int Miniball::nr_points () const
{
  return L.size();
}
   

Miniball::Cit Miniball::points_begin () const
{
  return L.begin();
}
   

Miniball::Cit Miniball::points_end () const
{
  return L.end();
} 

int Miniball::nr_support_points () const
{
  return B.support_size();
}
   

Miniball::Cit Miniball::support_points_begin () const
{
  return L.begin();
}
   

Miniball::Cit Miniball::support_points_end () const
{
  return support_end;
}

double Miniball::accuracy (double& slack) const
{
  double e, max_e = 0;
  int n_supp=0;
  Cit i;
  for (i=L.begin(); i!=support_end; ++i,++n_supp)
    if ((e = std::abs (B.excess (*i))) > max_e)
      max_e = e;
   
  // you've found a non-numerical problem if the following ever fails
  assert (n_supp == nr_support_points());
   
  for (i=support_end; i!=L.end(); ++i)
    if ((e = B.excess (*i)) > max_e)
      max_e = e;
   
  slack = B.slack();
  return (max_e/squared_radius());
}
   

bool Miniball::is_valid (double tolerance) const
{
  double slack;
  return ( (accuracy (slack) < tolerance) && (slack == 0) );
}   

// Miniball_b
// ----------
   

const double* Miniball_b::center () const
{
  return current_c;
}
   

double Miniball_b::squared_radius() const
{
  return current_sqr_r;
}
   

int Miniball_b::size() const
{
  return m;
}
   

int Miniball_b::support_size() const
{
  return s;
}
   

double Miniball_b::excess (const ANNpoint& p) const
{
  double e = -current_sqr_r;
  for (int k=0; k<d; ++k)
    e += mb_sqr(p[k]-current_c[k]);
  return e;
}

void Miniball_b::reset ()
{
  m = s = 0;
  // we misuse c[0] for the center of the empty sphere
  for (int j=0; j<d; ++j)
    c[0][j]=0;
  current_c = c[0];
  current_sqr_r = -1;
}

void Miniball_b::pop ()
{
  --m;
}

bool Miniball_b::push (const ANNpoint& p)
{
  int i, j;
  double eps = 1e-32;
  if (m==0) {
    for (i=0; i<d; ++i)
      q0[i] = p[i];
    for (i=0; i<d; ++i)
      c[0][i] = q0[i];
    sqr_r[0] = 0;
  } else {
    // set v_m to Q_m
    for (i=0; i<d; ++i)
      v[m][i] = p[i]-q0[i];
   
    // compute the a_{m,i}, i< m
    for (i=1; i<m; ++i) {
      a[m][i] = 0;
      for (j=0; j<d; ++j)
	a[m][i] += v[i][j] * v[m][j];
      a[m][i]*=(2/z[i]);
    }
   
    // update v_m to Q_m-\bar{Q}_m
    for (i=1; i<m; ++i) {
      for (j=0; j<d; ++j)
	v[m][j] -= a[m][i]*v[i][j];
    }
   
    // compute z_m
    z[m]=0;
    for (j=0; j<d; ++j)
      z[m] += mb_sqr(v[m][j]);
    z[m]*=2;
   
    // reject push if z_m too small
    if (z[m]<eps*current_sqr_r) {
      return false;
    }
   
    // update c, sqr_r
    double e = -sqr_r[m-1];
    for (i=0; i<d; ++i)
      e += mb_sqr(p[i]-c[m-1][i]);
    f[m]=e/z[m];
   
    for (i=0; i<d; ++i)
      c[m][i] = c[m-1][i]+f[m]*v[m][i];
    sqr_r[m] = sqr_r[m-1] + e*f[m]/2;
  }
  current_c = c[m];
  current_sqr_r = sqr_r[m];
  s = ++m;
  return true;
}
      

double Miniball_b::slack () const
{
  double* l = new double[d+1];
  double min_l=0;
  l[0] = 1;
  for (int i=s-1; i>0; --i) {
    l[i] = f[i];
    for (int k=s-1; k>i; --k)
      l[i]-=a[k][i]*l[k];
    if (l[i] < min_l) min_l = l[i];
    l[0] -= l[i];
  }
  if (l[0] < min_l) min_l = l[0];
  delete[] l;
  return ( (min_l < 0) ? -min_l : 0);
}
