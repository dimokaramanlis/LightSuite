C++/CUDA GPU-accelerated convolution in 2D and 3D. [![View GPU CUDA convolution 2D 3D on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/129964-gpu-cuda-convolution-2d-3d)

Based on NVIDIA cuda-samples convolutionFFT2D combined with matlab mexGPUexample.m.
- https://github.com/NVIDIA/cuda-samples/tree/master/Samples/5_Domain_Specific/convolutionFFT2D
- https://uk.mathworks.com/help/parallel-computing/run-mex-functions-containing-cuda-code.html 

I provide compiled .mexw64 files from a Windows 10 and compiled .mexa64 files from unix,
which should run out of the box.  
If this doesn't work for you due to different machine, a new mex compilation 
will be attempted and the
NVIDIA CUDA toolbox - including an nvcc compiler, supported C++ compiler, and library cuFFT - must be installed.

Run functions ``CUDAconvolution(data, kernel)`` or ``CUDAconvolution3D(data, kernel)`` analogous to matlab ``conv2``, ``convn``.

The method is convolution by FFT, pointwise multiply, and inverse FFT.

This method is much faster in the case of medium to large kernels; 
outperforms matlab starting at kernel size ~12 x 12 x 12 and speedup is more than 1000x at convolution 900x900x200 with 100x100x100 kernel (``test3d.mlx``).
Execution time should be constant and is <1s on my machine up to GPU memory limit.

Data should be the bigger array, as an array cut to original dimensions of data is returned.
