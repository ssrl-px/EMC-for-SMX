#include "orientCUDA.h"


#define gpuErr(ans) { gpuAssert((ans), __FILE__, __LINE__); }

__global__
//void orientMultiply(const double* __restrict__  qVecs, MAT3* rotMats, int numQ, int numRot, bool* out, double hcut, int minPred);
void orientMultiply(VEC3*  qVecs, MAT3* rotMats, int numQ, int numRot, bool* out, float hcut, int minPred);


inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void error_msg(cudaError_t err, const char* msg){
    if (err != cudaSuccess){
        printf("%s: CUDA error message: %s\n", msg, cudaGetErrorString(err));
        exit(err);
    }
}

void free_orientMatch(gpuOrient& gpu){
    if (gpu.rotMats != NULL)
        gpuErr(cudaFree(gpu.rotMats));
    if(gpu.qVecs != NULL)
        gpuErr(cudaFree(gpu.qVecs));
    if (gpu.out != NULL)
        gpuErr(cudaFree(gpu.out));
}

void setup_orientMatch(int dev_id, int maxNumQ, gpuOrient& gpu, np::ndarray Umats, bool alloc ){
    int numRot = Umats.shape(0)/9;
    if (alloc){
        gpu.numRot = numRot;
        gpu.max_numQ = maxNumQ;
        gpuErr(cudaSetDevice(dev_id));
        gpu.device = dev_id;
        gpuErr(cudaMallocManaged((void **)&gpu.rotMats, numRot*sizeof(MAT3)));
        gpuErr(cudaMallocManaged((void **)&gpu.out, numRot*sizeof(bool)));
        gpuErr(cudaMallocManaged((void **)&gpu.qVecs, maxNumQ*sizeof(VEC3)));
        //gpuErr(cudaMalloc((void **)&gpu.qVecs, maxNumQ*3*sizeof(double)));

        MAT3 Umat; // orientation matrix
        for (int i_rot=0; i_rot < numRot; i_rot ++){
            int i= i_rot*9;
            float uxx = float(bp::extract<double>(Umats[i]));
            float uxy = float(bp::extract<double>(Umats[i+1]));
            float uxz = float(bp::extract<double>(Umats[i+2]));
            float uyx = float(bp::extract<double>(Umats[i+3]));
            float uyy = float(bp::extract<double>(Umats[i+4]));
            float uyz = float(bp::extract<double>(Umats[i+5]));
            float uzx = float(bp::extract<double>(Umats[i+6]));
            float uzy = float(bp::extract<double>(Umats[i+7]));
            float uzz = float(bp::extract<double>(Umats[i+8]));
            Umat << uxx, uxy, uxz,
                    uyx, uyy, uyz,
                    uzx, uzy, uzz;
            gpu.rotMats[i_rot] = Umat;
        }
    }

    // optional size of each device block, else default to 128
    char* diffBragg_threads = getenv("ORIENT_THREADS_PER_BLOCK");
    if (diffBragg_threads==NULL)
        gpu.blockSize=128;
    else
        gpu.blockSize=atoi(diffBragg_threads);
    gpu.numBlocks = (numRot+gpu.blockSize-1)/gpu.blockSize;

}


void orientPeaks(gpuOrient& gpu, np::ndarray qvecs, float hcut, int minPred, bool verbose){

    double time;
    struct timeval t1, t2;//, t3 ,t4;

    gettimeofday(&t1, 0);
    int numQ = qvecs.shape(0)/3;

    if (verbose)
        printf("Setting cuda device %d\n", gpu.device);
    gpuErr(cudaSetDevice(gpu.device));
    if (numQ > gpu.max_numQ){
        printf("WARNING: re-allocating because maximum num Q vecs was exceeded!! Now maxNumQ =%d (was %d)\n",
            numQ, gpu.max_numQ);
        gpu.max_numQ = numQ;
        if (gpu.qVecs != NULL)
            gpuErr(cudaFree(gpu.qVecs));
        //gpuErr(cudaMalloc((void **)&gpu.qVecs, gpu.max_numQ*3*sizeof(double)));
        gpuErr(cudaMallocManaged((void **)&gpu.qVecs, gpu.max_numQ*sizeof(VEC3)));
    }

    // copy the Qvectors to the device
    //double* temp = new double[numQ*3];
    if (verbose)
        printf("Copying over %d qvectors to the GPU\n", numQ);
    for (int i_q=0; i_q < numQ; i_q++){
        int i = i_q*3;
        float qx = float(bp::extract<double>(qvecs[i]));
        float qy = float(bp::extract<double>(qvecs[i+1]));
        float qz = float(bp::extract<double>(qvecs[i+2]));
        //temp[i] = qx;
        //temp[i+1] = qy;
        //temp[i+2] = qz;
        VEC3 Q(qx,qy,qz);
        gpu.qVecs[i_q] = Q;
    }
    //gpuErr(cudaMemcpy(gpu.qVecs, temp, sizeof(double)*numQ*3, cudaMemcpyHostToDevice));
    //delete temp;
    gettimeofday(&t2, 0);
    time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;
    if(verbose)
        printf("Pre-kernel time=%f msec\n", time);

    gettimeofday(&t1, 0);
    // run the kernel
    orientMultiply<<<gpu.numBlocks, gpu.blockSize>>>(gpu.qVecs, gpu.rotMats, numQ, gpu.numRot, gpu.out, hcut, minPred);

    error_msg(cudaGetLastError(), "after kernel call");
    cudaDeviceSynchronize();
    gettimeofday(&t2, 0);
    time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;
    if(verbose)
        printf("kernel time=%f msec\n", time);

    gettimeofday(&t1, 0);
    bp::list rot_inds;
    for (int i=0; i< gpu.numRot; i++){
        if (gpu.out[i])
            rot_inds.append(i);
    }
    gpu.probable_rot_inds = rot_inds;
    gettimeofday(&t2, 0);
    time = (1000000.0*(t2.tv_sec-t1.tv_sec) + t2.tv_usec-t1.tv_usec)/1000.0;
    if(verbose)
        printf("POST kernel time=%f msec\n", time);

}


/* MAIN KERNEL */
__global__
void orientMultiply(VEC3* qVecs, MAT3* rotMats, int numQ, int numRot, bool* out, float hcut, int minPred){
//void orientMultiply(const double* __restrict__ qVecs, MAT3* rotMats, int numQ, int numRot, bool* out, double hcut, int minPred){

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int thread_stride = blockDim.x * gridDim.x;

    //extern __shared__ VEC3 s_qvecs[];
    //for(int i_q=0; i_q < numQ; i_q++){
    //    if (threadIdx.x==i_q){
    //        int i = i_q*3;
    //        double qx = __ldg(&qVecs[i]);
    //        double qy = __ldg(&qVecs[i+1]);
    //        double qz = __ldg(&qVecs[i+2]);
    //        VEC3 Q(qx,qy,qz);
    //        s_qvecs[i_q] =  Q;
    //    }
    //}
    //__syncthreads();

    for (int i_rot=tid; i_rot < numRot; i_rot += thread_stride){
        int count=0;
        for (int i_q=0; i_q < numQ; i_q ++ ){
            //int i=i_q*3;
            //double qx = __ldg(&qVecs[i]);
            //double qy = __ldg(&qVecs[i+1]);
            //double qz = __ldg(&qVecs[i+2]);
            //VEC3 Q(qx,qy,qz);
            //VEC3 Hkl = rotMats[i_rot]*Q;

            //VEC3 Hkl = rotMats[i_rot]*s_qvecs[i_q];
            VEC3 Hkl = rotMats[i_rot]*qVecs[i_q];

            //if (i_rot==0 && i_q==0){
            //    printf("h=%.3f,k=%.3f,l=%.3f\n", Hkl[0], Hkl[1], Hkl[2]);
            //    printf("qx=%.3f,qy=%.3f,qz=%.3f\n", Q[0], Q[1], Q[2]);
            //}
            //VEC3 Hkl = rotMats[i_rot]*qVecs[i_q];
            float h = ceil(Hkl[0]-0.5);
            float k = ceil(Hkl[1]-0.5);
            float l = ceil(Hkl[2]-0.5);
            VEC3 Hi(h,k,l);
            VEC3 deltaH = Hkl-Hi;
            float hnorm = deltaH.norm();
            if (hnorm < hcut)
                count += 1;
        }
        if (count >= minPred)
            out[i_rot] = true;
        else
            out[i_rot] = false;
    }
}



//cudaError_t CutlassSgemmNN(
//  int M,
//  int N,
//  int K,
//  float alpha,
//  float const *A,
//  int lda,
//  float const *B,
//  int ldb,
//  float beta,
//  float *C,
//  int ldc) {
//
//  // Define type definition for single-precision CUTLASS GEMM with column-major
//  // input matrices and 128x128x8 threadblock tile size (chosen by default).
//  //
//  // To keep the interface manageable, several helpers are defined for plausible compositions
//  // including the following example for single-precision GEMM. Typical values are used as
//  // default template arguments. See `cutlass/gemm/device/default_gemm_configuration.h` for more details.
//  //
//  // To view the full gemm device API interface, see `cutlass/gemm/device/gemm.h`
//
//  using ColumnMajor = cutlass::layout::ColumnMajor;
//
//  using CutlassGemm = cutlass::gemm::device::Gemm<float,        // Data-type of A matrix
//                                                  ColumnMajor,  // Layout of A matrix
//                                                  float,        // Data-type of B matrix
//                                                  ColumnMajor,  // Layout of B matrix
//                                                  float,        // Data-type of C matrix
//                                                  ColumnMajor>; // Layout of C matrix
//
//  // Define a CUTLASS GEMM type
//  CutlassGemm gemm_operator;
//
//  // Construct the CUTLASS GEMM arguments object.
//  //
//  // One of CUTLASS's design patterns is to define gemm argument objects that are constructible
//  // in host code and passed to kernels by value. These may include pointers, strides, scalars,
//  // and other arguments needed by Gemm and its components.
//  //
//  // The benefits of this pattern are (1.) a structured, composable strategy for passing host-constructible
//  // arguments to kernels and (2.) minimized initialization overhead on kernel entry.
//  //
//  CutlassGemm::Arguments args({M , N, K},  // Gemm Problem dimensions
//                              {A, lda},    // Tensor-ref for source matrix A
//                              {B, ldb},    // Tensor-ref for source matrix B
//                              {C, ldc},    // Tensor-ref for source matrix C
//                              {C, ldc},    // Tensor-ref for destination matrix D (may be different memory than source C matrix)
//                              {alpha, beta}); // Scalars used in the Epilogue
//
//  //
//  // Launch the CUTLASS GEMM kernel.
//  //
//
//  cutlass::Status status = gemm_operator(args);
//
//  //
//  // Return a cudaError_t if the CUTLASS GEMM operator returned an error code.
//  //
//
//  if (status != cutlass::Status::kSuccess) {
//    return cudaErrorUnknown;
//  }
//
//  // Return success, if no errors were encountered.
//  return cudaSuccess;
//}
//