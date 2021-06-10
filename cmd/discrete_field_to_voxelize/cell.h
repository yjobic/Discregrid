//
// Created by jobic on 09/06/2021.
//

#ifndef DISCREGRID_CELL_H
#define DISCREGRID_CELL_H

#include <vector>
#include <Discregrid/All>
#include <memory>

class cell {

  public:
    void setXdim(double,double);
    void setYdim(double,double);
    void setZdim(double,double);

    std::vector<double> getXdim();
    std::vector<double> getYdim();
    std::vector<double> getZdim();

  private:
    std::vector<double> xdim,ydim,zdim; //xmin,xmax,ymin,ymax,zmin,zmax defining the position of the cell
    bool border; //is this cell in a border solide/liquid ?
    bool solide; //is it a solid ?
    double sdfMeanValue; //mean value of the Signed distance function for this cell
};


#endif //DISCREGRID_CELL_H
