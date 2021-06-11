//
// Created by jobic on 10/06/2021.
//

template <typename T>
grid3D<T>::grid3D() {
  data=NULL;
  dimx=dimy=dimz=0;
}

template <typename T>
grid3D<T>::grid3D(const grid3D<T>& other) {
  unsigned int i,size=other.x*other.y*other.z;
  data=new T[size];
  dimx=other.w;
  dimy=other.h;
  dimz=other.l;
  for(i=0;i<size;++i){
    data[i]=other.data[i];
  }
}

template <typename T>
grid3D<T>::grid3D(unsigned int width,unsigned int height,unsigned int length) {
  dimx=width;
  dimy=height;
  dimz=length;
  data=new T[dimx*dimy*dimz];
}

template <typename T>
grid3D<T>::~grid3D(){
  delete[] data;
}

template <>
int grid3D<cell>::saveGrid(std::string filename) {
  if (filename.empty()) {
    std::cout << "empty filename" << std::endl;
    std::cout << "in saveGrid" << std::endl;
    exit(1);
  }

  std::ofstream outfile;
  outfile.open(filename.c_str());
  outfile << "VOXEL= " << std::endl;
  outfile << "NI= " << dimx << std::endl;
  outfile << "NJ= " << dimy << std::endl;
  outfile << "NK= " << dimz << std::endl;
  outfile << "END" << std::endl;
  outfile << "END" << std::endl;
  outfile << "ICON" << std::endl;
  for (unsigned int k = 0; k < dimz; ++k) {
    for (unsigned int j = 0; j < dimy; ++j) {
      for (unsigned int i = 0; i < dimx; ++i) {
        outfile << this->operator()(i,j,k).getSolide() << " ";
      }
      outfile << std::endl;
    }
    outfile << std::endl;
  }

  return 0;
}