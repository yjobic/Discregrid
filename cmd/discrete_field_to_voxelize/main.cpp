// run example options : -b "-5,5,-5,5,0,0.1" -r "20,20,20" -o out.toto ./small_cylindre.cdf
//
#include <Discregrid/All>
#include <Eigen/Dense>
#include <cxxopts/cxxopts.hpp>

#include <string>
#include <iostream>
#include <array>
#include "cell.h"
#include "grid3D.h"
#include "mesh.h"
#include "constants.h"

using namespace Eigen;


std::istream& operator>>(std::istream& is, std::array<unsigned int, 3>& data)
{
  is >> data[0] >> data[1] >> data[2];
  return is;
}

std::istream& operator>>(std::istream& is, std::array<double, 3>& data)
{
  is >> data[0] >> data[1] >> data[2] ;
  return is;
}

std::istream& operator>>(std::istream& is, std::array<double, 6>& data)
{
  is >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5] >> data[6];
  return is;
}

//resolution and sample is redondant, however i need it in order to check the sampling method of "integration"
// in order to determine if a cell is a fluid cell or a solid cell
void sdfToVTK(std::string filename, std::unique_ptr<Discregrid::DiscreteGrid> &sdf, unsigned int field_id,
              std::vector<unsigned int> resolution, std::vector<int> &sample, std::vector<double> bbox);


int main(int argc, char* argv[])
{
	cxxopts::Options options(argv[0], "Transforms a discrete SDF to a symetric voxel file.");
	options.positional_help("[input SDF file]");

	options.add_options()
	("h,help", "Prints this help text")
	("f,field_id", "ID in which the SDF to export is stored.", cxxopts::value<unsigned int>()->default_value("0"))
	("r,resolution", "Grid resolution", cxxopts::value<std::vector<unsigned int>>()->default_value("10,10,10"))
	("s,psymetry", "Symetry point", cxxopts::value<std::array<double, 3>>()->default_value("0. 0. 0."))
	("b,bbox", "bounding box", cxxopts::value<std::vector<double>>()->default_value("0.,0.,0.,0.,0.,0."))
	("n,nSample", "Number of samples per cell", cxxopts::value<std::vector<unsigned int>>()->default_value("10,10,10"))
	("o,output", "Output file", cxxopts::value<std::string>()->default_value(""))
	("input", "SDF file", cxxopts::value<std::vector<std::string>>())
	;

	try
	{
		options.parse_positional("input");
		auto result = options.parse(argc, argv);

		if (result.count("help"))
		{
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: SDFToBitmap -p xz file.sdf" << std::endl;
			exit(0);
		}
		if (!result.count("input"))
		{
			std::cout << "ERROR: No input file given." << std::endl;
			std::cout << options.help() << std::endl;
			std::cout << std::endl << std::endl << "Example: SDFToBitmap -p xz file.sdf" << std::endl;
			exit(1);
		}

		auto sdf = std::unique_ptr<Discregrid::DiscreteGrid>{};

		auto filename = result["input"].as<std::vector<std::string>>().front();
		auto lastindex = filename.find_last_of(".");
		auto extension = filename.substr(lastindex + 1, filename.length() - lastindex);

		std::cout << "Load SDF...";
		if (extension == "cdf" || extension == "cdm")
		{
			sdf = std::unique_ptr<Discregrid::CubicLagrangeDiscreteGrid>(
				new Discregrid::CubicLagrangeDiscreteGrid(filename));
		}
		std::cout << "DONE" << std::endl;

		auto const& domain = sdf->domain();
    auto resolution = result["r"].as<std::vector<unsigned int>>();
    auto psym = result["s"].as<std::array<double, 3>>();
		auto field_id = result["f"].as<unsigned int>();
    auto bbox = result["b"].as<std::vector<double>>();
    auto nbSamples = result["n"].as<std::vector<unsigned int>>();
    auto sample = std::vector<int>{(int)nbSamples[0],(int)nbSamples[1],(int)nbSamples[2]};

    if (!result.count("bbox")) {
      bbox= {domain.min()(0),domain.max()(0),
             domain.min()(1),domain.max()(1),
             domain.min()(2),domain.max()(2)};
    }

    //Creation of the mesh, considering the bounding box and the resolution
    mesh dm(bbox,resolution);


    dm.printResolution();
    dm.printBbox();

    std::cout << "symetry point : ";
    std::copy(std::begin(psym), std::end(psym), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;


    std::cout << "Nb samples per cell : ";
    std::copy(std::begin(sample), std::end(sample), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;


    // At first, we will not take into account the symetry point
    // Thus, the total number of cells known.
    grid3D<cell> grid (resolution[0],resolution[1],resolution[2]);

    std::cout << "Domaine total size : " << grid.size() << std::endl;
    dm.initCorrdsConstantStep();
    {
      std::vector<double> dims(6,0);

      for (unsigned int i=0; i < grid.getXdim(); ++i) {
        for (unsigned int j=0; j < grid.getYdim(); ++j) {
          for (unsigned int k=0; k < grid.getZdim(); ++k) {
            dims={dm.getXcoord(i),dm.getXcoord(i+1),
                  dm.getYcoord(j),dm.getYcoord(j+1),
                  dm.getZcoord(k),dm.getZcoord(k+1)};
            grid(i,j,k).setDims(dims);
          }
        }
      }
    }

    std::cout << "dims of (0,0,0) : " << std::endl;
	  grid(0,0,0).printDims();
    std::cout << std::endl;

    // init the fields of all the cells
    #pragma omp parallel for default(none) shared(resolution,field_id,sample,grid,sdf)
    for (unsigned int i=0; i < grid.getXdim(); ++i) {
      for (unsigned int j = 0; j < grid.getYdim(); ++j) {
        for (unsigned int k = 0; k < grid.getZdim(); ++k) {
          grid(i,j,k).initSdfMeanValue(sdf,field_id,sample);
        }
      }
    }

    std::cout << "Statistics of the first cell (0,0,0) " << std::endl;
    std::cout << "MeanValue : " << grid(0,0,0).getSdfMeanValue() << std::endl;
    std::cout << "ratioInsideOverTotal : " << grid(0,0,0).getRatioInsideOverTotal() << std::endl;

    std::cout << "Statistics of the first cell (1,1,0) " << std::endl;
    std::cout << "MeanValue : " << grid(1,1,0).getSdfMeanValue() << std::endl;
    std::cout << "ratioInsideOverTotal : " << grid(1,1,0).getRatioInsideOverTotal() << std::endl;

    auto ycellwidth = (bbox[YMAX]-bbox[YMIN])/resolution[1];
    auto zcellwidth = (bbox[ZMAX]-bbox[ZMIN])/resolution[2];
    double posX=0;
    double posZ=2*zcellwidth/sample[2];
    auto Pt = Vector3d{};
    Pt={posX,0.5-4*ycellwidth/9,posZ};
    std::cout << "value of point (" << Pt[0] << "," << Pt[1] << ","<< Pt[2] <<") : " << sdf->interpolate(field_id,Pt) << std::endl;
    Pt={posX,0.5-3*ycellwidth/9,posZ};
    std::cout << "value of point (" << Pt[0] << "," << Pt[1] << ","<< Pt[2] <<") : " << sdf->interpolate(field_id,Pt) << std::endl;
    Pt={posX,0.5-2*ycellwidth/9,posZ};
    std::cout << "value of point (" << Pt[0] << "," << Pt[1] << ","<< Pt[2] <<") : " << sdf->interpolate(field_id,Pt) << std::endl;
    Pt={posX,0.5-ycellwidth/9,posZ};
    std::cout << "value of point (" << Pt[0] << "," << Pt[1] << ","<< Pt[2] <<") : " << sdf->interpolate(field_id,Pt) << std::endl;
    Pt={posX,0.5,posZ};
    std::cout << "value of point (" << Pt[0] << "," << Pt[1] << ","<< Pt[2] <<") : " << sdf->interpolate(field_id,Pt) << std::endl;
    Pt={posX,0.5+ycellwidth/9,posZ};
    std::cout << "value of point (" << Pt[0] << "," << Pt[1] << ","<< Pt[2] <<") : " << sdf->interpolate(field_id,Pt) << std::endl;
    Pt={posX,0.5+2*ycellwidth/9,posZ};
    std::cout << "value of point (" << Pt[0] << "," << Pt[1] << ","<< Pt[2] <<") : " << sdf->interpolate(field_id,Pt) << std::endl;
    Pt={posX,0.5+3*ycellwidth/9,posZ};
    std::cout << "value of point (" << Pt[0] << "," << Pt[1] << ","<< Pt[2] <<") : " << sdf->interpolate(field_id,Pt) << std::endl;
    Pt={posX,0.5+4*ycellwidth/9,posZ};
    std::cout << "value of point (" << Pt[0] << "," << Pt[1] << ","<< Pt[2] <<") : " << sdf->interpolate(field_id,Pt) << std::endl;
    std::cout << std::endl;

    for (unsigned int i=0; i<sample[2]; ++i) {
      Pt={posX,0.5,bbox[ZMIN]+i*zcellwidth/((double)sample[2]-1)};
      std::cout << "value of point (" << Pt[0] << "," << Pt[1] << ","<< Pt[2] <<") : " << sdf->interpolate(field_id,Pt) << std::endl;
    }
    std::cout << std::endl;

    sdfToVTK("sdf.vtk", sdf, field_id, resolution, sample, bbox);
    grid.saveGrid("toto.txt");


	}
	catch (cxxopts::OptionException const& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
	
	return 0;
}

//resolution and sample is redondant, however i need it in order to check the sampling method of "integration"
// in order to determine if a cell is a fluid cell or a solid cell
void sdfToVTK(std::string filename, std::unique_ptr<Discregrid::DiscreteGrid> &sdf, unsigned int field_id,
              std::vector<unsigned int> resolution, std::vector<int> &sample, std::vector<double> bbox) {

  //save the SDF in a vtk file
  std::cout << "Save vtk for sdf" << std::endl;

  //seek the total number of point given the resolution and the sample size.
  //the resolution gives the number of cells, the sample size the number of samples in the cell.
  int sizeX =  resolution[0]*(sample[0]-1);
  int sizeY =  resolution[1]*(sample[1]-1);
  int sizeZ =  resolution[2]*(sample[2]-1);
  double xwidth=(bbox[XMAX]-bbox[XMIN])/sizeX;
  double ywidth=(bbox[YMAX]-bbox[YMIN])/sizeY;
  double zwidth=(bbox[ZMAX]-bbox[ZMIN])/sizeZ;
  sizeX++;sizeY++;sizeZ++;

  //we construct the points of the structure grid
  auto x=new double [sizeX];
  auto y=new double [sizeY];
  auto z=new double [sizeZ];
  for (unsigned int i=0; i<sizeX; ++i) {
    x[i] = bbox[XMIN] + i*xwidth;
  }
  for (unsigned int i=0; i<sizeY; ++i) {
    y[i] = bbox[YMIN] + i*ywidth;
  }
  for (unsigned int i=0; i<sizeZ; ++i) {
    z[i] = bbox[ZMIN] + i*zwidth;
  }

  auto evaluatedPoint = Vector3d{};
  grid3D<double> forVTK(sizeX,sizeY,sizeZ);
  #pragma omp parallel for default(none) shared(sizeX,sizeY,sizeZ,field_id,forVTK,sdf,x,y,z) private(evaluatedPoint)
  for (unsigned int k=0; k<sizeZ; ++k) {
    for (unsigned int j = 0; j < sizeY; ++j) {
      for (unsigned int i = 0; i < sizeX; ++i) {
        evaluatedPoint={x[i],y[j],z[k]};
        forVTK(i,j,k) = sdf->interpolate(field_id,evaluatedPoint);
      }
    }
  }

  std::cout << "sizeX " << sizeX << " sizeY " << sizeY << " sizeZ " << sizeZ << std::endl;
  std::cout << "xwidth " << xwidth << " ywidth " << ywidth << " zwidth " << zwidth << std::endl;

  //we write the file
  std::ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << std::endl;
  outfile << "Value of the signed distance function" << std::endl;
  outfile << "ASCII" << std::endl;
  outfile << "DATASET STRUCTURED_POINTS" << std::endl;
  outfile << "DIMENSIONS " << sizeX << " " << sizeY << " " << sizeZ << std::endl;
  outfile << "SPACING " << xwidth << " " << ywidth << " " << zwidth << std::endl;
  outfile << "ORIGIN " << bbox[XMIN] << " " << bbox[YMIN] << " " << bbox[ZMIN] << std::endl;
  outfile << "POINT_DATA " << sizeX*sizeY*sizeZ << std::endl;
  outfile << "SCALARS sdf double 1" << std::endl;
  outfile << "LOOKUP_TABLE default" << std::endl;
  for (unsigned int k=0; k<sizeZ; ++k) {
    for (unsigned int j = 0; j < sizeY; ++j) {
      for (unsigned int i = 0; i < sizeX; ++i) {
        outfile << forVTK(i,j,k) << " ";
      }
      outfile << std::endl;
    }
  }
  outfile.close();

  std::cout << "End of saving vtk" << std::endl;

}