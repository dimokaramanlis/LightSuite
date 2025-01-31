function dataOut = CUDAconvolution(data, kernel)
%CUDAconvolution(data, kernel) GPU-accelerated 2D convolution
% OS-specific mex file should be on path.
% Otherwise ,
% the code compiles mex file from mexGPUconvolution.cu, 
% convolutionFFT2D_common.h.  Then
%   - .cu and .h files must be on path
%   - NVIDIA nvcc compiler must be installed
%   - NVDIA cuFFT library must be installed
% INPUTS
% data, kernel - 2D numeric matlab arrays
% OUTPUTS
% dataOut - convolution of data and kernel, matlab 2D single array on host
%
% C++/CUDA code is adapted from NVIDIA cuda-samples/convolutionFFT2D
%--------------------------------------------------------------------------
%check inputs
assert(ismatrix(data) && ismatrix(kernel) && isnumeric(data) && isnumeric(kernel))
%--------------------------------------------------------------------------
% check compiled fun
ext = mexext;
if ~(exist(['mexGPUconvolution.' ext],"file")==3)
    disp("Compiling...");
    cufile = which('mexGPUconvolution.cu');
    eval(['mexcuda -v -lcufft ' cufile]);
end
%--------------------------------------------------------------------------
%check memory
if ~isgpuarray(data)
    GPUinfo = gpuDevice(1);
    GPUmemory= GPUinfo.AvailableMemory;
    sizeOfSingle =4;
    sizeOfComplexSingle=8;
    sizeOfData = numel(data);
    sizeOfPaddedData = prod(2.^ceil(log2(size(data)+size(kernel))));
    sizeOfKernel = numel(kernel);
    memoryNeeded =( (2*sizeOfData + sizeOfKernel + 2*sizeOfPaddedData) * sizeOfSingle + ...
        2*sizeOfPaddedData*sizeOfComplexSingle );
    if memoryNeeded >= GPUmemory
       error(['Inputs too large for GPU memory. \n' ...
           'Kernel size was %e elements, data size was %e elements.\n' ...
           'GPU memory to process would be %e bytes,' ...
           ' while GPU "%s" has %e bytes available.'], ...
           sizeOfKernel, sizeOfData, memoryNeeded, GPUinfo.Name, GPUmemory)
    end
    
    dataGPU=gpuArray(cast(data, 'single'));
    kernelGPU=gpuArray(cast(kernel, 'single'));
    [dataOutGPU] = mexGPUconvolution(dataGPU,kernelGPU);
    dataOut=gather(dataOutGPU);
else
    kernelGPU=gpuArray(cast(kernel, 'single'));
    [dataOutGPU] = mexGPUconvolution(data,kernelGPU);
    dataOut=gather(dataOutGPU);
end

end