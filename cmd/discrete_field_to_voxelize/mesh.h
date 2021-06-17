//
// Created by jobic on 16/06/2021.
//

#ifndef DISCREGRID_MESH_H
#define DISCREGRID_MESH_H

#include <vector>
#include <iostream>
#include <iterator>
#include <array>

class mesh {

  private:
    std::vector<double> Xcoords,Ycoords,Zcoords; // Coordinates of the voxels
    std::vector<unsigned int> resolution; // number of voxels in each directions {x,y,z}
    std::vector<double> bbox;    //bounding box of the mesh

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

    void printXcoords();
    void printYcoords();
    void printZcoords();
    void printResolution();
    void printBbox();
};


#endif //DISCREGRID_MESH_H
