
#include <Discregrid/All>
#include <Eigen/Dense>
#include <cxxopts/cxxopts.hpp>

#include <string>
#include <iostream>
#include <array>
#include "cell.h"
#include "grid3D.h"

using namespace Eigen;


std::istream& operator>>(std::istream& is, std::array<unsigned int, 3>& data)
{
  is >> data[0] >> data[1] >> data[2];
  return is;
}

std::istream& operator>>(std::istream& is, std::array<double, 3>& data)
{
  is >> data[0] >> data[1] >> data[2];
  return is;
}

int main(int argc, char* argv[])
{
	cxxopts::Options options(argv[0], "Transforms a discrete SDF to a symetric voxel file.");
	options.positional_help("[input SDF file]");

	options.add_options()
	("h,help", "Prints this help text")
	("f,field_id", "ID in which the SDF to export is stored.", cxxopts::value<unsigned int>()->default_value("0"))
	("r,resolution", "Grid resolution", cxxopts::value<std::array<unsigned int, 3>>()->default_value("10 10 10"))
	//("s,psymetry", "Symetry point", cxxopts::value<std::array<double, 3>>()->default_value("0. 0. 0."))
	//("s,psymetry", "Symetry point", cxxopts::value<std::array<double, 3>>()->default_value("0. 0. 0."))
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
    auto resolution = result["r"].as<std::array<unsigned int, 3>>();
    //auto psym = result["s"].as<std::vector<double>>();

		auto field_id = result["f"].as<unsigned int>();

		std::cout << "resolution : ";
    std::copy(std::begin(resolution), std::end(resolution), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;

    //std::cout << "symetry point : ";
    //std::copy(std::begin(psym), std::end(psym), std::ostream_iterator<double>(std::cout, " "));
    //std::cout << std::endl;

    // At first, we will not take into account the symetry point
    // Thus, the total number of cells known.
    grid3D<cell> grid (resolution[0],resolution[1],resolution[2]);

    std::cout << "Domaine total size : " << grid.size() << std::endl;
    {
      double xcellwidth,ycellwidth,zcellwidth,xmin,xmax,ymin,ymax,zmin,zmax;
      std::vector<double> dims(6,0);

      //simplification : the symetry point is at (0,0,0)
      //the widths are then
      xcellwidth = (domain.max()(0)-domain.min()(0))/resolution[0];
      ycellwidth = (domain.max()(0)-domain.min()(1))/resolution[1];
      zcellwidth = (domain.max()(0)-domain.min()(2))/resolution[2];

      std::cout << "xcellwidth : " << xcellwidth << std::endl;

      for (unsigned int i=0; i < resolution[0]; ++i) {
        for (unsigned int j=0; j < resolution[1]; ++j) {
          for (unsigned int k=0; k < resolution[2]; ++k) {
            xmin = domain.min()(0)+i*xcellwidth;
            xmax = domain.min()(0)+(i+1)*xcellwidth;
            ymin = domain.min()(1)+j*xcellwidth;
            ymax = domain.min()(1)+(j+1)*xcellwidth;
            zmin = domain.min()(2)+k*xcellwidth;
            zmax = domain.min()(2)+(k+1)*xcellwidth;

            dims={xmin,xmax,ymin,ymax,zmin,zmax};
            grid(i,j,k).setDims(dims);
          }
        }
      }
    }

    std::cout << "xdims of (0,0,0) : " << std::endl;
	  grid(0,0,0).printDims();
    std::cout << std::endl;

    auto sample = std::vector<int>{20,20,20};

    // init the fields of all the cells
    #pragma omp parallel for default(none) shared(resolution,field_id,sample,grid,sdf)
    for (unsigned int i=0; i < resolution[0]; ++i) {
      for (unsigned int j = 0; j < resolution[1]; ++j) {
        for (unsigned int k = 0; k < resolution[2]; ++k) {
          grid(i,j,k).initSdfMeanValue(sdf,field_id,sample);
        }
      }
    }

    //print the resulting boundary of the voxelized domain
    for (unsigned int k = 0; k < grid.getXdim(); ++k) {
      for (unsigned int j = 0; j < grid.getYdim(); ++j) {
        for (unsigned int i=0; i < grid.getZdim(); ++i) {
          std::cout << grid(i,j,k).getBorder() << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }

    grid.saveGrid("toto.txt");


	}
	catch (cxxopts::OptionException const& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
	
	return 0;
}