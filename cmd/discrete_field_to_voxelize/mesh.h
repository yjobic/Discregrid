//
// Created by jobic on 16/06/2021.
//

#ifndef DISCREGRID_MESH_H
#define DISCREGRID_MESH_H

#include <vector>
#include <iostream>
#include <iterator>

class mesh {

  public:
    mesh();
    ~mesh();
    void printXcoords();
    void printYcoords();
    void printZcoords();
    inline void addXcoords(double val) {Xcoords.push_back(val);}
    inline void addYcoords(double val) {Ycoords.push_back(val);}
    inline void addZcoords(double val) {Zcoords.push_back(val);}
    inline void addCoords(double valx, double valy, double valz) {
      Xcoords.push_back(valx);Ycoords.push_back(valy);Zcoords.push_back(valz);}

  private:
    std::vector<double> Xcoords,Ycoords,Zcoords;

};


#endif //DISCREGRID_MESH_H
