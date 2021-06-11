//
// Created by jobic on 09/06/2021.
//

#include "cell.h"

void cell::setDims(std::vector<double> &another) {
  assert(another.size()>0 && another.size()<=6);
  this->dims=another;
}

double cell::getDim(unsigned int i) {
  assert(i>0 && i<=6);
  return this->dims[i];
}

std::vector<double> cell::getDims() {
  return this->dims;
}
