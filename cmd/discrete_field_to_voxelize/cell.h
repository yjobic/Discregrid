//
// Created by jobic on 09/06/2021.
//

#ifndef DISCREGRID_CELL_H
#define DISCREGRID_CELL_H

#include <vector>
#include <Discregrid/All>
#include <memory>
#include "constants.h"


class cell {

  public:
    cell();
    ~cell();
    void setDims(std::vector<double> &another);
    void initSdfMeanValue(std::unique_ptr<Discregrid::DiscreteGrid> &sdf, unsigned int fields_id, std::vector<int> &samples);

    double getDim(unsigned int i);
    std::vector<double> getDims();
    inline int getSolide(){return solide;}
    inline int getBorder(){return border;}
    inline double getSdfMeanValue(){return sdfMeanValue;}
    inline double getRatioInsideOverTotal(){return ratioInsideOverTotal;}
    void printDims();

  private:
    std::vector<double> dims; //xmin,xmax,ymin,ymax,zmin,zmax defining the position of the cell
    int border; //is this cell in a border solide/liquid ?
    int solide; //is it a solid ?
    double sdfMeanValue; //mean value of the Signed distance function for this cell
    double ratioInsideOverTotal; //SDF is not very accurate for openBoundary domain. We know it's quite good
                                 //for closed one, which is the inside. Thus, computing the ratio of inside samples over
                                 //the total gives a better indicator of Solide/Fluide voxel.
};


#endif //DISCREGRID_CELL_H
