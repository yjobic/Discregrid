//
// Created by jobic on 09/06/2021.
//

#include "cell.h"
#include <Eigen/Dense>

using namespace Eigen;

cell::cell() {
  border=-1;
  solide=-1;
  sdfMeanValue=-123456;
}

cell::~cell() {
  border=-2;
  solide=-2;
  sdfMeanValue=-123456;
}

void cell::setDims(std::vector<double> &another) {
  assert(another.size()==6);
  this->dims=another;
}

double cell::getDim(unsigned int i) {
  assert(i>0 && i<=6);
  return this->dims[i];
}

std::vector<double> cell::getDims() {
  return this->dims;
}

void cell::printDims() {
  assert(dims.size()==6);
  std::copy(std::begin(dims), std::end(dims), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
}

void cell::initSdfMeanValue(std::unique_ptr<Discregrid::DiscreteGrid> &sdf, unsigned int fields_id, std::vector<int> &samples) {
  double  xwidths,ywidths,zwidths,dist,moyDist;
  bool initn=false,initp=false;
  auto sample = Vector3d{};

  assert(samples.size()==3);
  assert(dims.size()==6);

  xwidths = (dims[XMAX]-dims[XMIN])/samples[0];
  ywidths = (dims[YMAX]-dims[YMIN])/samples[1];
  zwidths = (dims[ZMAX]-dims[ZMIN])/samples[2];

  moyDist=0;
  for(unsigned int i=0; i<samples[0]; ++i) {
    for(unsigned int j=0; j<samples[1]; ++j) {
      for(unsigned int k=0; k<samples[2]; ++k) {
        sample(0) = dims[XMIN]+i*xwidths;
        sample(1) = dims[YMIN]+j*ywidths;
        sample(2) = dims[ZMIN]+k*zwidths;
        dist = sdf->interpolate(fields_id,sample);
        if (!initp && dist > 0) {
          initp=true;
        }
        if (!initn && dist < 0) {
          initn=true;
        }
        moyDist += dist;
      }
    }
  }
  sdfMeanValue = moyDist/(samples[0]+samples[1]+samples[2]);
  if (initp && initn) border=1;
  else border=0;
  if (sdfMeanValue>1e-13) solide=1;
  else solide=0;
}