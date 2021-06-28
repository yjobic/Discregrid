//
// Created by jobic on 16/06/2021.
//

#ifndef DISCREGRID_MESH_H
#define DISCREGRID_MESH_H

#include <vector>
#include <iostream>
#include <iterator>
#include <array>
#include "constants.h"
#include "cell.h"
#include "grid3D.h"


class mesh {

  private:
    std::vector<double> Xcoords,Ycoords,Zcoords; // Coordinates of the voxels
    std::vector<unsigned int> resolution; // number of voxels in each directions {x,y,z}
    std::vector<double> bbox;    //bounding box of the mesh
    grid3D<cell> grid;

  public:
    mesh(std::vector<double> bbox, std::vector<unsigned int> resolution);
    ~mesh();
    inline void addXcoords(double val) {Xcoords.push_back(val);}
    inline void addYcoords(double val) {Ycoords.push_back(val);}
    inline void addZcoords(double val) {Zcoords.push_back(val);}
    inline void addCoords(double valx, double valy, double valz) {
      Xcoords.push_back(valx);Ycoords.push_back(valy);Zcoords.push_back(valz);}

    inline double getXcoord(unsigned int pos) {return Xcoords[pos];}
    inline double getYcoord(unsigned int pos) {return Ycoords[pos];}
    inline double getZcoord(unsigned int pos) {return Zcoords[pos];}

    void initCorrdsConstantStep();

    // Getters / Setters of the grid
    inline void printGridSize() {grid.size();}

    void printXcoords();
    void printYcoords();
    void printZcoords();
    void printResolution();
    void printBbox();
};


#endif //DISCREGRID_MESH_H
