//
// Created by jobic on 10/06/2021.
//

#ifndef DISCREGRID_GRID3D_H
#define DISCREGRID_GRID3D_H

#include <iostream>

template<typename T>
class grid3D {

  protected:
    T* data;
    unsigned int dimx,dimy,dimz;

  public:
    grid3D(){
      data=NULL;
      dimx=dimy=dimz=0;}
    grid3D(const grid3D<T>& other){
      unsigned int i,size=other.x*other.y*other.z;
      data=new T[size];
      dimx=other.w;
      dimy=other.h;
      dimz=other.l;
      for(i=0;i<size;++i){data[i]=other.data[i];}}
    grid3D(unsigned int width,unsigned int height,unsigned int length){
      dimx=width;
      dimy=height;
      dimz=length;
      data=new T[dimx*dimy*dimz];}
    ~grid3D(){
      delete[] data;}
    inline T& operator()(unsigned int x,unsigned int y,unsigned int z){
      return data[z + dimz*(y+dimy*x)];}
    inline const T& operator()(unsigned int x,unsigned int y,unsigned int z) const{
      return data[z + dimz*(y+dimy*x)];}
    inline unsigned int size() const{
      return dimx*dimy*dimz;}
};


#endif //DISCREGRID_GRID3D_H
