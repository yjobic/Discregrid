//
// Created by jobic on 09/06/2021.
//

#include "cell.h"

void cell::setXdim(double min, double max) {
  this->xdim.push_back(min);
  this->xdim.push_back(max);
}

void cell::setYdim(double min, double max) {
  this->xdim.push_back(min);
  this->xdim.push_back(max);
}

void cell::setZdim(double min, double max) {
  this->xdim.push_back(min);
  this->xdim.push_back(max);
}

std::vector<double> cell::getXdim() {
  return this->xdim;
}

std::vector<double> cell::getYdim() {
  return this->ydim;
}

std::vector<double> cell::getZdim() {
  return this->zdim;
}

