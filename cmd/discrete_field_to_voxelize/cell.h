//
// Created by jobic on 09/06/2021.
//

#ifndef DISCREGRID_CELL_H
#define DISCREGRID_CELL_H

#include <vector>
#include <Discregrid/All>
#include <memory>

#define XMIN 0
#define XMAX 1
#define YMIN 2
#define YMAX 3
#define ZMIN 4
#define ZMAX 5

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
    void printDims();

  private:
    std::vector<double> dims; //xmin,xmax,ymin,ymax,zmin,zmax defining the position of the cell
    int border; //is this cell in a border solide/liquid ?
    int solide; //is it a solid ?
    double sdfMeanValue; //mean value of the Signed distance function for this cell
};


#endif //DISCREGRID_CELL_H
