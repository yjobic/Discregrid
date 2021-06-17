//
// Created by jobic on 16/06/2021.
//

#include "mesh.h"

mesh::mesh(std::vector<double> bboxInput, std::vector<unsigned int> resolutionInput) {
  if (bboxInput.size() != 6) {
    std::cout << "bounding box must be of dim 6!" << std::endl;
    exit(1);
  }

  if (resolutionInput.size() != 3) {
    std::cout << "Number of voxels (resolution) must be of dim 3!" << std::endl;
    exit(1);
  }

  bbox=bboxInput;
  resolution=resolutionInput;
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

void mesh::printResolution() {
  std::cout << "resolution : ";
  std::copy(std::begin(resolution), std::end(resolution), std::ostream_iterator<unsigned int>(std::cout, " "));
  std::cout << std::endl;
}

void mesh::printBbox() {
  std::cout << "bounding box : ";
  std::copy(std::begin(bbox), std::end(bbox), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
}