/*
 * Example of how to use the mxGPUArray API in a MEX file.  This example shows
 * how to write a MEX function that takes a gpuArray input and returns a
 * gpuArray output, e.g. dataOut=mexFunction(data).
 *
 * Copyright 2012 The MathWorks, Inc.
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"

// Include CUDA runtime and CUFFT
#include <cuda_runtime.h>
#include <cufft.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "convolutionFFT2D_common.h"

#define USE_TEXTURE 1
#define POWER_OF_TWO 1

#if (USE_TEXTURE)
#define LOAD_FLOAT(i) tex1Dfetch<float>(texFloat, i)
#define SET_FLOAT_BASE
#else
#define LOAD_FLOAT(i) d_Src[i]
#define SET_FLOAT_BASE
#endif

////////////////////////////////////////////////////////////////////////////////
/// Position convolution kernel center at (0, 0) in the image
////////////////////////////////////////////////////////////////////////////////
__global__ void padKernel_kernel(float *d_Dst, float *d_Src, int fftH, int fftW,
                                 int kernelH, int kernelW, int kernelY,
                                 int kernelX
#if (USE_TEXTURE)
                                 ,
                                 cudaTextureObject_t texFloat
#endif
                                 ) {

  const int y = blockDim.y * blockIdx.y + threadIdx.y;
  const int x = blockDim.x * blockIdx.x + threadIdx.x;

  if ( y < kernelH && x < kernelW) {

    int ky = y - kernelY;

    if (ky < 0) {
      ky += fftH;
    }

    int kx = x - kernelX;

    if (kx < 0) {
      kx += fftW;
    }

    d_Dst[ ky * fftW + kx] = LOAD_FLOAT(y * kernelW + x);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Prepare data for "pad to border" addressing mode
////////////////////////////////////////////////////////////////////////////////
__global__ void unpad_kernel(float *d_Dst, float *d_Src,
                                             int fftH, int fftW, int dataH,
                                            int dataW
#if (USE_TEXTURE)
                                            ,
                                            cudaTextureObject_t texFloat
#endif
                                            ) {
  const int y = blockDim.y * blockIdx.y + threadIdx.y;
  const int x = blockDim.x * blockIdx.x + threadIdx.x;

  if (y < dataH && x < dataW) {

    d_Dst[y * dataW + x] = LOAD_FLOAT(y * fftW + x);
  }
}

__global__ void padDataClampToBorder_kernel(float *d_Dst, float *d_Src,
                                             int fftH, int fftW, int dataH,
                                            int dataW,  int kernelH, int kernelW,
                                            int kernelY, int kernelX
#if (USE_TEXTURE)
                                            ,
                                            cudaTextureObject_t texFloat
#endif
                                            ) {
  const int y = blockDim.y * blockIdx.y + threadIdx.y;
  const int x = blockDim.x * blockIdx.x + threadIdx.x;
  const int borderH = dataH + kernelY;
  const int borderW = dataW + kernelX;

  if (y < fftH && x < fftW) {
    int  dy, dx;

    if (y < dataH) {
      dy = y;
    }

    if (x < dataW) {
      dx = x;
    }

    if (y >= dataH && y < borderH) {
      dy = dataH - 1;
    }

    if (x >= dataW && x < borderW) {
      dx = dataW - 1;
    }

    if (y >= borderH) {
      dy = 0;
    }

    if (x >= borderW) {
      dx = 0;
    }

    d_Dst[y * fftW + x] = LOAD_FLOAT(dy * dataW + dx);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Modulate Fourier image of padded data by Fourier image of padded kernel
// and normalize by FFT size
////////////////////////////////////////////////////////////////////////////////
inline __device__ void mulAndScale(fComplex &a, const fComplex &b,
                                   const float &c) {
  fComplex t = {c * (a.x * b.x - a.y * b.y), c * (a.y * b.x + a.x * b.y)};
  a = t;
}

__global__ void modulateAndNormalize_kernel(fComplex *d_Dst, fComplex *d_Src,
                                            int dataSize, float c) {
  const int i = blockDim.x * blockIdx.x + threadIdx.x;

  if (i >= dataSize) {
    return;
  }

  fComplex a = d_Src[i];
  fComplex b = d_Dst[i];

  mulAndScale(a, b, c);

  d_Dst[i] = a;
}



////////////////////////////////////////////////////////////////////////////////
/// Position convolution kernel center at (0, 0) in the image
////////////////////////////////////////////////////////////////////////////////
extern "C" void padKernel(float *d_Dst, float *d_Src, int fftH, int fftW,
                           int kernelH, int kernelW, int kernelY, int kernelX) {
  assert(d_Src != d_Dst);
  dim3 threads( 8, 4);
  dim3 grid(iDivUp(kernelW, threads.x), iDivUp(kernelH, threads.y));

  SET_FLOAT_BASE;
#if (USE_TEXTURE)
  cudaTextureObject_t texFloat;
  cudaResourceDesc texRes;
  memset(&texRes, 0, sizeof(cudaResourceDesc));

  texRes.resType = cudaResourceTypeLinear;
  texRes.res.linear.devPtr = d_Src;
  texRes.res.linear.sizeInBytes = sizeof(float) * kernelH * kernelW ;
  texRes.res.linear.desc = cudaCreateChannelDesc<float>();

  cudaTextureDesc texDescr;
  memset(&texDescr, 0, sizeof(cudaTextureDesc));

  texDescr.normalizedCoords = false;
  texDescr.filterMode = cudaFilterModeLinear;
  texDescr.addressMode[0] = cudaAddressModeWrap;
  texDescr.readMode = cudaReadModeElementType;

  cudaCreateTextureObject(&texFloat, &texRes, &texDescr, NULL);
#endif

  padKernel_kernel<<<grid, threads>>>(d_Dst, d_Src, fftH, fftW, 
                                      kernelH, kernelW, kernelY, kernelX
#if (USE_TEXTURE)
                                      ,
                                      texFloat
#endif
                                      );

#if (USE_TEXTURE)
  cudaDestroyTextureObject(texFloat);
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Prepare data for "pad to border" addressing mode
////////////////////////////////////////////////////////////////////////////////
extern "C" void padDataClampToBorder(float *d_Dst, float *d_Src, 
                                     int fftH, int fftW, int dataH, int dataW,
                                     int kernelH, int kernelW,
                                     int kernelY, int kernelX) {
  assert(d_Src != d_Dst);
  dim3 threads(8, 8);
  dim3 grid(iDivUp(fftW, threads.x), iDivUp(fftH, threads.y));

#if (USE_TEXTURE)
  cudaTextureObject_t texFloat;
  cudaResourceDesc texRes;
  memset(&texRes, 0, sizeof(cudaResourceDesc));

  texRes.resType = cudaResourceTypeLinear;
  texRes.res.linear.devPtr = d_Src;
  texRes.res.linear.sizeInBytes = sizeof(float) * dataH * dataW ;
  texRes.res.linear.desc = cudaCreateChannelDesc<float>();

  cudaTextureDesc texDescr;
  memset(&texDescr, 0, sizeof(cudaTextureDesc));

  texDescr.normalizedCoords = false;
  texDescr.filterMode = cudaFilterModeLinear;
  texDescr.addressMode[0] = cudaAddressModeWrap;
  texDescr.readMode = cudaReadModeElementType;

  cudaCreateTextureObject(&texFloat, &texRes, &texDescr, NULL);
#endif

  padDataClampToBorder_kernel<<<grid, threads>>>(
      d_Dst, d_Src,  fftH, fftW,  dataH, dataW,  kernelH, kernelW, kernelY, kernelX
#if (USE_TEXTURE)
      ,
      texFloat
#endif
      );

#if (USE_TEXTURE)
  cudaDestroyTextureObject(texFloat);
#endif
}


extern "C" void unpad(float *d_Dst, float *d_Src, 
                                     int fftH, int fftW, int dataH, int dataW ) {
  assert(d_Src != d_Dst);
  dim3 threads(8, 8);
  dim3 grid(iDivUp(fftW, threads.x), iDivUp(fftH, threads.y));

#if (USE_TEXTURE)
  cudaTextureObject_t texFloat;
  cudaResourceDesc texRes;
  memset(&texRes, 0, sizeof(cudaResourceDesc));

  texRes.resType = cudaResourceTypeLinear;
  texRes.res.linear.devPtr = d_Src;
  texRes.res.linear.sizeInBytes = sizeof(float) * fftH * fftW ;
  texRes.res.linear.desc = cudaCreateChannelDesc<float>();

  cudaTextureDesc texDescr;
  memset(&texDescr, 0, sizeof(cudaTextureDesc));

  texDescr.normalizedCoords = false;
  texDescr.filterMode = cudaFilterModeLinear;
  texDescr.addressMode[0] = cudaAddressModeWrap;
  texDescr.readMode = cudaReadModeElementType;

  cudaCreateTextureObject(&texFloat, &texRes, &texDescr, NULL);
#endif

  unpad_kernel<<<grid, threads>>>(
      d_Dst, d_Src,  fftH, fftW,  dataH, dataW
#if (USE_TEXTURE)
      ,
      texFloat
#endif
      );

#if (USE_TEXTURE)
  cudaDestroyTextureObject(texFloat);
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Modulate Fourier image of padded data by Fourier image of padded kernel
// and normalize by FFT size
////////////////////////////////////////////////////////////////////////////////
extern "C" void modulateAndNormalize(fComplex *d_Dst, fComplex *d_Src, 
                                     int fftH, int fftW, int padding) {
  assert(fftW % 2 == 0);
  const int dataSize =  fftH *(fftW / 2 + padding);
  
  modulateAndNormalize_kernel<<<iDivUp(dataSize, 256), 256 >>>(
      d_Dst, d_Src, dataSize, 1.0f / (float)(fftW * fftH));
}



/*
 * Host code
 */


int snapTransformSize(int dataSize) {
  int hiBit;
  unsigned int lowPOT, hiPOT;

  dataSize = iAlignUp(dataSize, 16);

  for (hiBit = 31; hiBit >= 0; hiBit--)
    if (dataSize & (1U << hiBit)) {
      break;
    }

  lowPOT = 1U << hiBit;

  if (lowPOT == (unsigned int)dataSize) {
    return dataSize;
  }

  hiPOT = 1U << (hiBit + 1);

  if (hiPOT <= 1024) {
    return hiPOT;
  } else {
    return iAlignUp(dataSize, 512);
  }
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    /* Declare all variables.*/
    mxGPUArray const *data;
    mxGPUArray const *kernel;
    //detect dimensions of data

    float *d_Data;
    float *d_Kernel;
    mxGPUArray *PaddedData_c;
    mxGPUArray *PaddedKernel_c;
    float *d_PaddedKernel;
    float *d_PaddedData;
    mxGPUArray *DataSpectrum_c, * KernelSpectrum_c;
    fComplex* d_DataSpectrum, * d_KernelSpectrum;
    cufftHandle fftPlanFwd, fftPlanInv;
    float *d_dataOut;
    mxGPUArray *dataOut_c;
    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.  Must be matlab single type.";
    char const * const errMsg3D = "Invalid input to MEX file.  3D arrays expected.";

    /* Choose a reasonably sized number of threads for the block. */
    int const threadsPerBlock = 256;
    int blocksPerGrid;

    /* Initialize the MathWorks GPU API. */
    mxInitGPU();

    /* Throw an error if the input is not a GPU array. */
    if ((nrhs!=2) || !(mxIsGPUArray(prhs[0]))|| !(mxIsGPUArray(prhs[1]))) {
        mexErrMsgIdAndTxt(errId, errMsg);
    }

    data = mxGPUCreateFromMxArray(prhs[0]);
    kernel = mxGPUCreateFromMxArray(prhs[1]);

    //check inputs
    if ((mxGPUGetClassID(data) != mxSINGLE_CLASS)||(mxGPUGetClassID(kernel) != mxSINGLE_CLASS)) { //goes with float
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    if ((mxGPUGetNumberOfDimensions(data) != 2)||(mxGPUGetNumberOfDimensions(kernel) != 2) ) { 
        mexErrMsgIdAndTxt(errId, errMsg3D);
    }

    
    //detect input dimensions
    mwSize const * const dimsData = mxGPUGetDimensions(data);
    mwSize const dataW = dimsData[0]; //I want W to be innermost dimension
    mwSize const dataH = dimsData[1];
    mwSize const * const dimsKernel = mxGPUGetDimensions(kernel);
    mwSize const kernelW = dimsKernel[0]; 
    mwSize const kernelH = dimsKernel[1];
    mwSize const kernelX = kernelW/2; 
    mwSize const kernelY = kernelH/2;
    unsigned int const fftW = snapTransformSize(dataW + kernelW - 1);
    unsigned int const fftH = snapTransformSize(dataH + kernelH - 1);
    mwSize const dimsFft[2] = {fftW, fftH};
    mwSize const dimsComplex[2] = { fftW/2+1 , fftH};

    /*
     * Now that we have verified the data type, extract a pointer to the input
     * data on the device.
     */
    d_Data = (float *)(mxGPUGetDataReadOnly(data));
    d_Kernel = (float *)(mxGPUGetDataReadOnly(kernel));
    /* Create GPUArray on device only. */
    PaddedData_c = mxGPUCreateGPUArray(2,
                            dimsFft,
                            mxGPUGetClassID(data),
                            mxGPUGetComplexity(data),
                            MX_GPU_INITIALIZE_VALUES);
    d_PaddedData = (float *)(mxGPUGetData(PaddedData_c));
    PaddedKernel_c = mxGPUCreateGPUArray(2,
                            dimsFft,
                            mxGPUGetClassID(kernel),
                            mxGPUGetComplexity(kernel),
                            MX_GPU_INITIALIZE_VALUES);
    d_PaddedKernel = (float *)(mxGPUGetData(PaddedKernel_c));
    //fourier space complex arrays
    KernelSpectrum_c = mxGPUCreateGPUArray(2,
                            dimsComplex,
                            mxGPUGetClassID(data), 
                            mxCOMPLEX, 
                            MX_GPU_INITIALIZE_VALUES);
    DataSpectrum_c = mxGPUCreateGPUArray(2,
                            dimsComplex,
                            mxGPUGetClassID(data),
                            mxCOMPLEX,
                            MX_GPU_INITIALIZE_VALUES);
    d_KernelSpectrum = (fComplex *)(mxGPUGetData(KernelSpectrum_c));
    d_DataSpectrum = (fComplex *)(mxGPUGetData(DataSpectrum_c));
    dataOut_c = mxGPUCreateGPUArray(2,
                            dimsData,
                            mxGPUGetClassID(data),
                            mxGPUGetComplexity(data),
                            MX_GPU_DO_NOT_INITIALIZE);
    d_dataOut = (float *)(mxGPUGetData(dataOut_c));
    


    //pad data and kernel
    padKernel(d_PaddedKernel, d_Kernel, fftH, fftW, kernelH, kernelW, kernelY,
        kernelX);
    padDataClampToBorder(d_PaddedData, d_Data, fftH, fftW, dataH, dataW,
        kernelH, kernelW, kernelY, kernelX);

    //Fourier transform
    cufftPlan2d(&fftPlanFwd, fftH, fftW, CUFFT_R2C);
    cufftPlan2d(&fftPlanInv, fftH, fftW, CUFFT_C2R);
    cufftExecR2C(fftPlanFwd, (cufftReal*)d_PaddedKernel,
        (cufftComplex*)d_KernelSpectrum);
    cufftExecR2C(fftPlanFwd, (cufftReal*)d_PaddedData,
        (cufftComplex*)d_DataSpectrum);

    // multiply elementwise in fourier space + normalize
    modulateAndNormalize((fComplex*)d_DataSpectrum, (fComplex*)d_KernelSpectrum,  fftH, fftW, 1);

    //inverse Fourier transform
    cufftExecC2R(fftPlanInv, (cufftComplex*)d_DataSpectrum,
        (cufftReal*)d_PaddedData);

    //unpad result

    unpad(d_dataOut, d_PaddedData, fftH, fftW,dataH, dataW);

    /* Wrap the result up as a MATLAB gpuArray for return. */
//     plhs[0] = mxGPUCreateMxArrayOnGPU(kernel);
//     plhs[1] = mxGPUCreateMxArrayOnGPU(data);
//     mxGPUArray * fftkernel = mxGPUCopyReal(KernelSpectrum_c);
//     plhs[2] = mxGPUCreateMxArrayOnGPU(fftkernel);
//     mxGPUArray * fftkernel_i = mxGPUCopyImag(KernelSpectrum_c);
//     plhs[3] = mxGPUCreateMxArrayOnGPU(fftkernel_i);
//     mxGPUArray * fftconvolved = mxGPUCopyReal(DataSpectrum_c);
//     plhs[4] = mxGPUCreateMxArrayOnGPU(fftconvolved);
//     mxGPUArray * fftconvolved_i = mxGPUCopyImag(DataSpectrum_c);
//     plhs[5] = mxGPUCreateMxArrayOnGPU(fftconvolved_i);
//     plhs[6] = mxGPUCreateMxArrayOnGPU(PaddedKernel_c);
    plhs[0] = mxGPUCreateMxArrayOnGPU(dataOut_c);

    /*
     * The mxGPUArray pointers are host-side structures that refer to device
     * data. These must be destroyed before leaving the MEX function.
     */
    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);

    mxGPUDestroyGPUArray(data);
    mxGPUDestroyGPUArray(kernel);
    mxGPUDestroyGPUArray(PaddedKernel_c);
    mxGPUDestroyGPUArray(PaddedData_c);
    mxGPUDestroyGPUArray(KernelSpectrum_c);
    mxGPUDestroyGPUArray(DataSpectrum_c);
}