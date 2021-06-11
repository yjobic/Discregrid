
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

    auto xsamples = 5; auto ysamples = xsamples; auto zsamples = xsamples;

    auto data = std::vector<double>{};
		data.resize(xsamples * ysamples * zsamples);

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
      std::vector<double> dims;

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
	  auto vecToPrintStdin = grid(0,0,0).getDims();
    std::copy(std::begin(vecToPrintStdin), std::end(vecToPrintStdin), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    //std::cout << "xdims of (1,2,3) : " << grid(1,2,3).getXdim()[0] << grid(1,2,3).getXdim()[1] << std::endl;
    //std::cout << "xdims of (2,3,4) : " << grid(2,3,4).getXdim()[0] << grid(2,3,4).getXdim()[1] << std::endl;

//#pragma omp parallel for
//		for (int k = 0; k < static_cast<int>(xsamples * ysamples); ++k)
//		{
//			auto i = k % xsamples;
//			auto j = k / xsamples;
//
//			auto xr = static_cast<double>(i) / static_cast<double>(xsamples);
//			auto yr = static_cast<double>(j) / static_cast<double>(ysamples);
//
//			auto x = domain.min()(dir(0)) + xr * diag(dir(0)) + 0.5 * xwidth;
//			auto y = domain.min()(dir(1)) + yr * diag(dir(1)) + 0.5 * ywidth;
//
//			auto sample = Vector3d{};
//			sample(dir(0)) = x;
//			sample(dir(1)) = y;
//			sample(dir(2)) = domain.min()(dir(2)) + 0.5 * (1.0 + depth) * diag(dir(2));
//
//			data[k] = sdf->interpolate(field_id, sample);
//			if (data[k] == std::numeric_limits<double>::max())
//			{
//				data[k] = 0.0;
//			}
//		}
//
//		std::cout << "DONE" << std::endl;
//
//		auto min_v = *std::min_element(data.begin(), data.end());
//		auto max_v = *std::max_element(data.begin(), data.end());
//
//		auto out_file = result["o"].as<std::string>();
//		if (out_file == "")
//		{
//			out_file = filename;
//			if (out_file.find(".") != std::string::npos)
//			{
//				auto lastindex = out_file.find_last_of(".");
//				out_file = out_file.substr(0, lastindex);
//			}
//			out_file += ".bmp";
//		}
//
//		std::cout << "Ouput file: " << out_file << std::endl;
//
//		std::cout << "Export BMP...";
//		std::transform(data.begin(), data.end(), data.begin(), [&max_v, &min_v](double v) {return v >= 0.0 ? v / std::abs(max_v) : v / std::abs(min_v); });
//
//		auto pixels = std::vector<std::array<unsigned char, 3u>>(data.size());
//
//		auto cm = result["c"].as<std::string>();
//		if (cm != "gb" && cm != "rs")
//		{
//			std::cerr << "WARNING: Unknown color map option. Fallback to mode 'gb'." << std::endl;
//		}
//
//		if (cm == "gb")
//			std::transform(data.begin(), data.end(), pixels.begin(), doubleToGreenBlueInverse);
//		else if (cm == "rs")
//			std::transform(data.begin(), data.end(), pixels.begin(), doubleToRedSequential);
//
//		BmpReaderWriter::saveFile(out_file.c_str(), xsamples, ysamples, &pixels.front()[0]);
//		std::cout << "DONE" << std::endl;
//
//		std::cout << std::endl << "Statistics:" << std::endl;
//		std::cout << "\tdomain         = " << domain.min().transpose() << ", " << domain.max().transpose() << std::endl;
//		std::cout << "\tmin value      = " << min_v << std::endl;
//		std::cout << "\tmax value      = " << max_v << std::endl;
//		std::cout << "\tbmp resolution = " << xsamples << " x " << ysamples << std::endl;
	}
	catch (cxxopts::OptionException const& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
	
	return 0;
}