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
  grid=grid3D<cell>(resolutionInput[0],resolutionInput[1],resolutionInput[2]);
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

void mesh::initCorrdsConstantStep() {
  double xcellwidth,ycellwidth,zcellwidth,xmin,xmax,ymin,ymax,zmin,zmax;

  //simplification : the symetry point is at (0,0,0)
  //the widths are then
  xcellwidth = (bbox[XMAX]-bbox[XMIN])/resolution[0];
  ycellwidth = (bbox[YMAX]-bbox[YMIN])/resolution[1];
  zcellwidth = (bbox[ZMAX]-bbox[ZMIN])/resolution[2];

  std::cout << "xcellwidth : " << xcellwidth <<
            " ycellwidth : " << ycellwidth <<
            " zcellwidth : " << zcellwidth << std::endl;

  for (unsigned int i=0; i < resolution[0]; ++i) {
    xmin = bbox[XMIN] + i * xcellwidth;
    xmax = bbox[XMIN] + (i + 1) * xcellwidth;
    addXcoords(xmin);
    if (i == resolution[0]-1) addXcoords(xmax);
  }

  for (unsigned int j=0; j < resolution[1]; ++j) {
    ymin = bbox[YMIN] + j * ycellwidth;
    ymax = bbox[YMIN] + (j + 1) * ycellwidth;
    addYcoords(ymin);
    if (j == resolution[1]-1) addYcoords(ymax);
  }

  for (unsigned int k=0; k < resolution[2]; ++k) {
    zmin = bbox[ZMIN] + k * zcellwidth;
    zmax = bbox[ZMIN] + (k + 1) * zcellwidth;
    addZcoords(zmin);
    if (k == resolution[2] - 1) addZcoords(zmax);
  }
}