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
  ratioInsideOverTotal=-1;
}

cell::~cell() {
  border=-2;
  solide=-2;
  sdfMeanValue=-123456;
  ratioInsideOverTotal=-1;
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
  int     nbInside;
  double  xwidths,ywidths,zwidths,dist,moyDist;
  bool initn=false,initp=false;
  auto sample = Vector3d{};

  assert(samples.size()==3);
  assert(dims.size()==6);

  xwidths = (dims[XMAX]-dims[XMIN])/samples[0];
  ywidths = (dims[YMAX]-dims[YMIN])/samples[1];
  zwidths = (dims[ZMAX]-dims[ZMIN])/samples[2];

  nbInside=0;
  for(unsigned int i=0; i<samples[0]; ++i) {
    for(unsigned int j=0; j<samples[1]; ++j) {
      for(unsigned int k=0; k<samples[2]; ++k) {
        sample(0) = dims[XMIN]+i*xwidths;
        sample(1) = dims[YMIN]+j*ywidths;
        sample(2) = dims[ZMIN]+k*zwidths;
        dist = sdf->interpolate(fields_id,sample);
        if (dist == std::numeric_limits<double>::max())  {
          dist = 0.0;
        }
        if (!initp && dist > 0) {
          initp=true;
        }
        if (!initn && dist < 0) {
          initn=true;
        }
        moyDist += dist;
        if (dist>=0) nbInside++;
      }
    }
  }
  unsigned int totalSamples = samples[0]*samples[1]*samples[2];
  sdfMeanValue = moyDist/totalSamples;
  //SDF is not very accurate for openBoundary domain. We know it's quite good
  //for closed one, which is the inside. Thus, computing the ratio of inside samples over
  //the total gives a better indicator of Solide/Fluide voxel.
  ratioInsideOverTotal = nbInside/(double)totalSamples;
  if (ratioInsideOverTotal>=0.5) solide=1;
  else solide=0;
//  if (moyDist>0) solide=1;
//  else solide=0;
  if (initp && initn) border=1;
  else border=0;
}