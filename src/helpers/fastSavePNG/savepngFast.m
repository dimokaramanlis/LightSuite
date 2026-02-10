
function savepngFast(hfig, path, name, varargin)
%
%%% savepngFast %%%
%
%
% This function save the figure to png file much faster than MATLAB -print
% routine with good quality.
%
% ===============================Inputs====================================
%
%   hfig : handle to the figure (default: gcf).
%   path : path for the target .png file (default: pwd).
%   name : name of the target .png file (default: figure1).
%   resolution : quality of .png file (default: 600 DPI).
%   magnify : relevant to quality of output (default : 2).
%   renderer : renderering method (default: -opengl).
%   compression : A number between 0 and 10 controlling the amount of
%                compression to try to achieve with PNG file. 0 implies no
%                compresson, fastest option. 10 implies the most amount of
%                compression, slowest option (default: 7).
%
%================================Output====================================
%
%   .png file
%
%   Note : this function uses print2array function from export_fig function
%   group and some other internal function from same package to read the
%   image into data array. It also uses savepng function and mex file to
%   save the figure.
%
% written by Mohammad, 14.09.2015
%
%
if nargin < 1
    hfig = gcf;
end;

if nargin < 2
    path = pwd;
end;

if nargin < 3
    name = 'figure1';
end

if nargin > 3
    resolution = varargin{1};
else
    resolution = 600;
end;

if nargin > 4
    magnify = varargin{2};
else
    magnify = resolution/300;
end

if nargin > 5
    renderer = varargin{3};
else
    renderer = 1;
end

if nargin > 6
    compression = varargin{4};
else
    compression = 7;
end

% Set the renderer
switch renderer
    case 1
        renderer = '-opengl';
    case 2
        renderer = '-zbuffer';
    case 3
        renderer = '-painters';
    otherwise
        renderer = '-opengl'; % Default for bitmaps
end

% getting the image data
cData = print2array(hfig, magnify, renderer);

% resolution = ['-r',num2str(resolution)];
% magnify = ['-m',num2str(magnify)];
% cData = export_fig(hfig,resolution,magnify,renderer,'-a1','-nocrop');

fullname = [path,'\',name,'.png'];
% saving to png
savepng(cData,fullname, compression, resolution);
