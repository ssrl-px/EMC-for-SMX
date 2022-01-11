
#include <Eigen/Dense>
typedef Eigen::Matrix<float,3,1> VEC3;
typedef Eigen::Matrix<float,3,3> MAT3;
#include <iostream>
#include <boost/python/numpy.hpp>
#include <stdio.h>
#include <sys/time.h>

//#include <cutlass/cutlass.h>
//#include <cutlass/numeric_types.h>
//#include <cutlass/core_io.h>
//#include <cutlass/gemm/device/gemm.h>
//#include <cutlass/util/host_tensor.h>


namespace bp = boost::python;
namespace np = boost::python::numpy;

struct gpuOrient {
  MAT3* rotMats=NULL;
  VEC3* qVecs=NULL;
  //double* qVecs=NULL;
  bool* out=NULL;
  int max_numQ;
  int numRot;
  int numBlocks, blockSize;
  int device;
  bp::list probable_rot_inds;
};

void orientPeaks(gpuOrient& gpu,
                np::ndarray qvecs,
                float hcut,
                int minpred, bool verbose);
void setup_orientMatch(int dev_id, int maxNumQ, gpuOrient& gpu, np::ndarray Umats, bool alloc);
void free_orientMatch(gpuOrient& gpu);
