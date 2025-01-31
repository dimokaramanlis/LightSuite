function blobs = blobdetect3D(image, diameter,varargin)
%BLOBDETECT3D Find locations of dark/bright features.
%

%argparse
p = inputParser;

addRequired(p,'image',@(x) isnumeric(x)&&ndims(x)==3);
addRequired(p,'diameter', @(x) isempty(x)||(isnumeric(x)&&0<x));
addParameter(p,'Filter', 'LoG', @(x)ischar(x));
addParameter(p,'DarkBackground', true, @(x)islogical(x));
addParameter(p,'MedianFilter', true, @(x)islogical(x));
addParameter(p,'KernelSize', [], @(x)isnumeric(x)&&x>=1);
addParameter(p,'OverlapFilter', false, @(x)islogical(x));
addParameter(p,'QualityFilter', 0.1, @(x)isnumeric(x)&&x>=0 &&x<=1);
addParameter(p,'BorderWidth', [], @(x)isnumeric(x)&&x>=0); 
addParameter(p,'GPU', false, @(x)islogical(x));

parse(p, image, diameter, varargin{:})

if isempty(p.Results.BorderWidth)
   BorderWidth=ceil(diameter/2);
else
   BorderWidth=p.Results.BorderWidth;
end

if p.Results.DarkBackground
    for i=1:size(image,3)
        %TODO there has to be a better way than loop
        image(:,:,i)=imcomplement(image(:,:,i));
    end
end

if p.Results.MedianFilter
    %with matlab default connectivity and radius
    image=medfilt3(image);
end
%make the LoG kernel
sigma_ = diameter/(sqrt(3)*2); %sigma = r / sqrt(D) 

if isempty(p.Results.KernelSize)
    %taken from trackmate DetectionUtils.java createLoGKernel ,
    %creditted Tobias Gauss
    halfKernelSize = max(2,ceil(3*sigma_+.5));
    kernelSize = 3+2*halfKernelSize;
else
    %if user input is even, will be effectively rounded up to
    %nearest odd in use
    kernelSize = p.Results.KernelSize;
end

%check on GPU
if p.Results.GPU
    if gpuDeviceCount==0
        warning("No GPU devices found.  Running as 'GPU=false'");
        p.Results.GPU =false;
    end
end

kernel=LoG_kernel_3D(sigma_, kernelSize);

%convolve with the kernel
if p.Results.GPU
    try
        image_conv = CUDAconvolution3D(image, kernel);
    catch e
        disp(e.message);
        disp("CUDAconvolution3D not found or not working, " + ...
            "Falling back to matlab pieced convolution");
        try
            image_conv = convn(gpuArray(image), gpuArray(kernel), 'same');
            image_conv = gather(image_conv);
        catch
            image_conv = convn(image, kernel, 'same');
        end
    end
else
    try
        image_conv = convn(gpuArray(image), gpuArray(kernel), 'same');
        image_conv = gather(image_conv);
    catch
        image_conv = convn(image, kernel, 'same');
    end
end

%detect local maxima
%alternative taken from scipy: apply maximum filter (=imdilate).
% take points where original image == filtered image to be maxima
SE=strel('cube',4); %squares should be more separable -faster?
maxima = image_conv==imdilate(image_conv,SE);

%handle edges, erase maxima near edge
maxima([1:BorderWidth end-BorderWidth:end],:,:)=0;
maxima(:,[1:BorderWidth end-BorderWidth:end],:)=0;
try
    maxima(:,:,[1:BorderWidth end-BorderWidth:end])=0;
catch Error
    disp('Maybe you entered a 2D image into blobdetect3D')
    rethrow(Error)
end

%return list of coordinates
[xs,ys, zs]=find3D(maxima);

% each maximum is assiged a quality score.  for now just intensity.
quality = zeros([size(xs),1]);
for i=1:size(xs)
    quality(i)=image_conv(xs(i), ys(i), zs(i));
end
quality=rescale(quality);

if ~isempty(p.Results.QualityFilter)
    q_min=p.Results.QualityFilter;
    xs=xs(quality>q_min);
    ys=ys(quality>q_min);
    zs=zs(quality>q_min);
    quality=quality(quality>q_min);
end
%3D version of overlap filter
if p.Results.OverlapFilter
    [xs,ys, zs, quality] = overlap_filter3D(xs,ys,zs,quality,diameter);
end

    blobs = [xs,ys,zs, quality];
    if ~isempty(blobs)
    blobs = array2table(blobs,...
    VariableNames={'x','y','z','Intensity'});
    else
        %return empty table
        blobs= table('Size',[0,4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'x','y','z','Intensity'});
    end
end

function [xs, ys, zs] = find3D(array)
[xs, y_s] = find(array);
[ys, zs] = ind2sub([size(array,2), size(array, 3)], y_s);
end

function kernel = Laplacian_kernel_3D()
%simple kernel for testing
kernel=zeros([3 3 3]);
kernel(:,:,1)= [0,0,0;0,-1,0;0,0,0];
kernel(:,:,2)= [0,-1,0;-1,6,-1;0,-1,0];
kernel(:,:,3)= [0,0,0;0,-1,0;0,0,0];
end   

function kernel = LoG_kernel_3D(sigma_, kernel_size)
%arrays of x, y coordinates distance from center
%handle odd/even kernel_size - always make an odd kernel size,
%rounding up
range=-floor(kernel_size/2):floor(kernel_size/2);
[x,y,z] = ndgrid(range,range, range);
%LoG kernel formula
d=(x.^2+y.^2+z.^2)/(2*sigma_^2);
kernel= -1/((pi*2)^(3/2)*sigma_^5)*(3-2*d).*exp(-d);
%TODO should I discretize to ints?
end

function kernel = Gaussian_kernel_nD(sigma_sq, kernel_size)
%DoG faster?  Good for small objects.
end

function [xs, ys,zs, quality] = overlap_filter3D(xs,ys,zs,quality,diameter)
del_list = [];
rad_sq=(diameter/2)^2;
for i=1:size(xs)
    if ~ismember(i, del_list)
    for j=i+1:size(xs)
        if ~ismember(j, del_list) && (xs(i)-xs(j))^2+(ys(i)-ys(j))^2+(zs(i)-zs(j))^2 < rad_sq
            %delete weaker object - add to list to remove later
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
zs(del_list)=[];
quality(del_list)=[];
end
