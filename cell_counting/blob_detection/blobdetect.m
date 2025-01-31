function blobs = blobdetect(image, diameter,varargin)
%BLOBDETECT Find locations of dark/bright features.
%

%argparse
p = inputParser;

addRequired(p,'image',@(x) isnumeric(x)&&ismatrix(x));
addRequired(p,'diameter', @(x) isempty(x)||(isnumeric(x)&&0<x));

addParameter(p,'Filter', "LoG", @(x) any(validatestring(x,{'LoG', 'DoG'})) );
addParameter(p,'DarkBackground', true, @(x)islogical(x));
addParameter(p,'MedianFilter', true, @(x)islogical(x));
addParameter(p,'OverlapFilter', false, @(x)islogical(x));
addParameter(p,'QualityFilter', 0.1, @(x)isnumeric(x)&&x>=0 &&x<=1);
addParameter(p,'KernelSize', [], @(x)isnumeric(x)&&x>=1);
addParameter(p,'IntensityFilter', true, @(x)islogical(x));
addParameter(p,'BorderWidth', [], @(x)isnumeric(x)&&x>=0);
addParameter(p,'GPU', false, @(x)islogical(x));

parse(p, image, diameter, varargin{:})

if p.Results.DarkBackground
    image=imcomplement(image);
end

if p.Results.MedianFilter
    %using matlab's default connectivity and radius is
    image=medfilt2(image);
end

if isempty(p.Results.BorderWidth)
   BorderWidth=ceil(diameter/2);
else
   BorderWidth=p.Results.BorderWidth;
end

switch p.Results.Filter
    case "LoG"
        %make the LoG kernel
        sigma = (diameter/(sqrt(2)^3)); %sigma = r / sqrt(D)
        if isempty(p.Results.KernelSize)
            %taken from trackmate DetectionUtils.java createLoGKernel ,
            %creditted Tobias Gauss
            halfKernelSize = max(2,ceil(3*sqrt(sigma)+.5));
            kernelSize = 3+2*halfKernelSize;
        else
            %if user input is even, will be effectively rounded up to
            %nearest odd in use
            kernelSize = p.Results.KernelSize;
        end
        kernel=LoG_kernel(sigma, kernelSize);
        %convolve with the kernel
        if p.Results.GPU
            try
                image_conv = CUDAconvolution(image, kernel);
            catch e
                disp(e.message);
                disp("CUDAconvolution not found or not working, " + ...
                    "Falling back to matlab pieced convolution");
                try
                    image_conv = conv2(gpuArray(image), gpuArray(kernel), 'same');
                    image_conv = gather(image_conv);
                catch
                    image_conv = conv2(image, kernel, 'same');
                end
            end
        else
            try
                image_conv = conv2(gpuArray(image), gpuArray(kernel), 'same');
                image_conv = gather(image_conv);
            catch
                image_conv = convn(image, kernel, 'same');
            end
        end
    case "DoG"
        %smaller DoG sigma^2
        sigma_sq_1 = (diameter/(sqrt(2)^3)*0.9)^2;
        %larger DoG sigma^2
        sigma_sq_2 = (diameter/(sqrt(2)^3)*1.1)^2;
        kernel_1=Gaussian_kernel(sigma_sq_1, kernelSize);
        kernel_2=Gaussian_kernel(sigma_sq_2, kernelSize);
        %FFTconvolve with each kernel
        image_conv_1 = conv2(image, kernel_1, 'same');
        image_conv_2 = conv2(image, kernel_2, 'same');
        image_conv = image_conv_2 - image_conv_1 ;
end
% detect local maxima
% method taken from scipy: apply maximum filter (=imdilate). 
% take points where original image == filtered image to be maxima
SE=strel('square',4); 
maxima = image_conv==imdilate(image_conv,SE);

%handle edge
maxima(1:BorderWidth,:)=0;
maxima(end-BorderWidth:end,:)=0;
maxima(:,end-BorderWidth:end)=0;
maxima(:,1:BorderWidth)=0;

%return list of coordinates
[xs,ys]=find(maxima);

%TODO each maximum is assiged a quality score.  for now just intensity.
quality = zeros([size(xs, 1),1]);
for i=1:size(xs, 1)
    quality(i)=image_conv(xs(i),ys(i));
end
quality=rescale(quality);

if ~isempty(p.Results.QualityFilter)
    q_min=p.Results.QualityFilter;
    xs=xs(quality>q_min);
    ys=ys(quality>q_min);
    quality=quality(quality>q_min);
end

if p.Results.OverlapFilter
    [xs,ys, quality] = overlap_filter(xs,ys,quality,diameter);
end
% return output
blobs = [ys,xs,quality];

end

function kernel = Laplacian_kernel()
%simple kernel for testing
kernel= [0,-1,0;-1,4,-1;0,-1,0];
end

function kernel = LoG_kernel(sigma, kernel_size)
%arrays of x, y coordinates distance from center
%handle odd/even kernel_size - always makes it odd kernel size,
%rounding up
range=-floor(kernel_size/2):floor(kernel_size/2);
[x,y] = ndgrid(range,range);
%LoG kernel formula
d=(x.^2+y.^2)/(2*sigma^2);
kernel= -1/(pi*sigma^4)*(1-d).*exp(-d);
end

function kernel = Gaussian_kernel(sigma_sq, kernel_size)
%for DoG
range=-floor(kernel_size/2):floor(kernel_size/2);
[x,y] = ndgrid(range,range);
%Just a gaussian kernel (2D)
d=(x.^2+y.^2)/(2*sigma_sq);
kernel= 1/(2*pi*sigma_sq)*exp(-d);
end

function [xs, ys, quality] = overlap_filter(xs,ys,quality,diameter)
del_list = []; %list of indices to delete later
rad_sq=(diameter/2)^2;
for i=1:size(xs)
    if ~ismember(i, del_list)
        for j=i+1:size(xs)
            if ~ismember(j, del_list) && (xs(i)-xs(j))^2+(ys(i)-ys(j))^2 < rad_sq
                %delete weaker object
                if quality(i)<quality(j)
                    del_list(end+1)=i;
                else
                    del_list(end+1)=j;
                end
            end
        end
    end
end
xs(del_list)=[];
ys(del_list)=[];
quality(del_list)=[];
end
