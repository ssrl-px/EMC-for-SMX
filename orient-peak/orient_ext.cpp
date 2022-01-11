#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <iostream>
#include "orientCUDA.h"
#define BOOST_LIB_NAME "boost_numpy"
#include <boost/config/auto_link.hpp>

namespace bp=boost::python;
namespace np=boost::python::numpy;

class probaOr{
public:
  virtual ~probaOr(){}
  // constructor
  probaOr(){}
  gpuOrient gpu;
  inline void alloc(int device_id, np::ndarray rotations, int maxQvecs){
    int num_rot=rotations.shape(0)/9;
    printf("Determined number of rotations=%d\n", num_rot);
    //printf("We will allocate space for %d orientations!\n", num_rot);
    setup_orientMatch( device_id, maxQvecs, gpu, rotations, true);
  }
  inline void free(){
    free_orientMatch(gpu);
  }
  inline void oriPeaks(np::ndarray qvecs,
                       float hcut, int minWithinHcut, bool verbose){
    orientPeaks(gpu, qvecs, hcut, minWithinHcut, verbose);
  }

  inline bp::list listOrients(){
    return gpu.probable_rot_inds;
  }

  inline void print_rotMat(int i_rot){
    MAT3 M = gpu.rotMats[i_rot];
    printf("Rotation matrix %d=\n%.7f %.7f %.7f\n%.7f %.7f %.7f\n%.7f %.7f %.7f\n",
       M(0,0), M(0,1), M(0,2),
       M(1,0), M(1,1), M(1,2),
       M(2,0), M(2,1), M(2,2));
  }

};

BOOST_PYTHON_MODULE(orient){
    Py_Initialize();
    np::initialize();
    typedef bp::return_value_policy<bp::return_by_value> rbv;
    typedef bp::default_call_policies dcp;
    typedef bp::return_internal_reference<> rir;
    bp::class_<probaOr>("probable_orients", bp::no_init)
      //.def(bp::init<double>((bp::arg_("val")),"ORIENTS"))
      .def(bp::init<>("returns a class instance"))
      .def ("allocate_orientations", &probaOr::alloc, "move the orientations to the device")
      .def ("orient_peaks", &probaOr::oriPeaks, "compute probable orientations (main CUDA kernel)")
      .def("free_device", &probaOr::free, "free any allocated GPU memory")
      .def ("print_rotMat", &probaOr::print_rotMat, "show elements of allocated rotMat i_rot")
      .def ("get_probable_orients", &probaOr::listOrients, "returns a list of rotation matrix indices")
      //.add_property("num_ori",
      //               make_getter(&probaOr::num_orient,rbv()),
      //               make_setter(&probaOr::num_orient,dcp()),
      //               "Number of orientations.")
    ;
}