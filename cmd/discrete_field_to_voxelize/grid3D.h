//
// Created by jobic on 10/06/2021.
//

#ifndef DISCREGRID_GRID3D_H
#define DISCREGRID_GRID3D_H

#include <iostream>
#include <fstream>
#include <string>

template <typename T>
class grid3D {

  protected:
    T* data;
    unsigned int dimx,dimy,dimz;

  public:
    grid3D();
    grid3D(const grid3D<T>& other);
    grid3D(unsigned int width,unsigned int height,unsigned int length);
    ~grid3D();
    inline unsigned int getXdim(){return dimx;}
    inline unsigned int getYdim(){return dimy;}
    inline unsigned int getZdim(){return dimz;}
    int saveGrid(std::string filename);

    inline T& operator()(unsigned int x,unsigned int y,unsigned int z){
      return data[z + dimz*(y+dimy*x)];}
    inline const T& operator()(unsigned int x,unsigned int y,unsigned int z) const{
      return data[z + dimz*(y+dimy*x)];}
    inline unsigned int size() const{
      return dimx*dimy*dimz;}
};

#include "grid3D.tpp"
#endif //DISCREGRID_GRID3D_H
