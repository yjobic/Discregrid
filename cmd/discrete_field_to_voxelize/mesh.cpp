//
// Created by jobic on 16/06/2021.
//

#include "mesh.h"

mesh::mesh() {
  Xcoords=Ycoords=Zcoords={};
}

mesh::~mesh() {

}

void mesh::printXcoords() {
  std::cout << "Xcoords" << std::endl;
  std::copy(std::begin(Xcoords), std::end(Xcoords), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  std::cout << std::endl;
}

void mesh::printYcoords() {
  std::cout << "Ycoords" << std::endl;
  std::copy(std::begin(Ycoords), std::end(Ycoords), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  std::cout << std::endl;
}

void mesh::printZcoords() {
  std::cout << "Zcoords" << std::endl;
  std::copy(std::begin(Zcoords), std::end(Zcoords), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  std::cout << std::endl;
}