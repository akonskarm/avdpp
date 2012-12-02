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

//SCALE all point inside e/11.0 hypersphere
//TODO: try to remove the stupid vector
//#include "globals.h"
#include "scale.h"
//#include <iostream>
//#include <fstream>
//#include "Miniball_dynamic_d.h"
void scale()
{
  Miniball mb(Globals::dim);
  std::string buffer;
  //std::ifstream infile ("./random_numbers");
  std::ifstream infile ("./real_data");

  for (int i=0; i!=Globals::len; i++) { 
    Globals::points.push_back(annAllocPt(Globals::dim));
    for (int j=0; j!=Globals::dim; j++) {
      getline(infile,buffer);
      Globals::points[i][j] = atof(buffer.c_str());
    }
    
    //TODO:comment these 3 lines to use random_numbers
    for (int j=0; j!=7-Globals::dim; j++) {
      getline(infile,buffer);
    }

    mb.check_in(Globals::points[i]);    
  }
  mb.build();
  
//  Globals::scale = Globals::epsilon / (sqrt(mb.squared_radius()) * 11.0);
//  std::cout << "SCALE " << Globals::scale << std::endl;
  Globals::scale = pow(2,floor(log2(Globals::epsilon / (sqrt(mb.squared_radius()) * 11.0))));

//  std::cout << "SCALE " << Globals::scale << std::endl;
  
  Globals::move = annAllocPt(Globals::dim);
  ANNpoint temp = mb.center();
  for (int i=0; i!=Globals::dim; i++) {
    //Globals::move[i] = fabs(temp[i] * Globals::scale - 0.5);
    Globals::move[i] = temp[i] * Globals::scale - 0.5;
  }
//  std::cout << "MOVE ";
//  printPt(std::cout, Globals::move);
  
  for (int i=0; i<Globals::len; i++) {
    for (int j=0; j<Globals::dim; j++) {
      Globals::points[i][j] *= Globals::scale;
      Globals::points[i][j] -= Globals::move[j];
      //printPt(std::cout, Globals::points[i]);
      //std::cout << Globals::points[i][j] << std::endl;
      assert(Globals::points[i][j]*0.95 <= 0.5 + Globals::epsilon / 11.0);
      assert(Globals::points[i][j]*1.05 >= 0.5 - Globals::epsilon / 11.0);
    }
  }
  delete[] temp;
  infile.close();
  /*for (int i=0; i<Globals::len; i++)
  {
    printPt(std::cout,Globals::points[i]);
  }*/
  //delete temp??
  //delete mb???
}
