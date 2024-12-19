function [newvol,opts] = prepareRegistrationImage(opts,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    if isfield(opts, 'regvolpath')
        backvol = load(opts.regvolpath);
    else
        dirtocheck = fileparts(opts.fproc);
        candfile   = dir(fullfile(dirtocheck, 'temporary*.mat'));
        backvol    = load(fullfile(candfile.folder, candfile.name));
    end
    backvol    = backvol.voldown;
else
    backvol = varargin{1};
end
% for determing edges (mask)
% regpath   = fullfile(opts.savepath, 'registration_volume.tif');

%------------------------------------------------------------------
% fix uniformity of background
% for background imadjust(imflatfield(backsignal, 1000))
% findchangepts(std(backvol(:,:,500),[],1),'MaxNumChanges',2,'Statistic','rms')
% save background

% sigmause  = size(backvol, 2)/10;
% isample   = round(linspace(size(backvol,3)*0.3, size(backvol,3)*0.7, 50));
% yr        = round(size(backvol, 1)/4);
% xr        = round(size(backvol, 2)/3);
% yuse      = round(size(backvol, 1)/2)+(-yr:yr);
% xuse      = round(size(backvol, 2)/2)+(-xr:xr);
% 
% for ii = 1:numel(isample)
%     backvol(:, :, isample(ii))
%     profcurr = medfilt2(backvol(yuse, xuse, isample(1)), [121 121]);
% end
% 
%%
% figure out a mask
% Noutx     = round(size(backvol, 2) * 0.1);
% Nouty     = round(size(backvol, 1) * 0.1);
% quantbase = ceil(quantile(backvol([1:Nouty end-Nouty:end], [1:Noutx end-Noutx:end], :), 0.99, 'all'));
% maskuse   = backvol > quantbase;
% maskuse   = medfilt3(maskuse, 5*[1 1 1]);
% 

%%
centpx    = round(size(backvol)/2);
naround   = round(min(centpx)/3);
centind   = round(size(backvol)/2) + (-naround:naround)';

% normalize contrast left-right and turn to uint8
% basevol   = imresize3(backvol, 0.5);
% basevol   = gpuArray(basevol);
% tic;
% bg = imgaussfilt3(basevol, 25, 'FilterDomain','spatial');
% toc;
% basevol = basevol./bg;
% basevol = imresize3(gather(basevol), size(backvol));

bvalinds  = randperm(numel(backvol), 1e5);

topval    = quantile(backvol(centind(:,1), centind(:,2),centind(:,3)), 0.999,'all')*2;
bottomval = quantile(backvol(bvalinds), 0.01, 'all');
newvol    = (backvol-bottomval)/(topval-bottomval);
newvol    = uint8(newvol * 255);
%--------------------------------------------------------------------------
% save as tif volume
opts.regvolsize = size(newvol);
makeNewDir(opts.savepath)
save(fullfile(opts.savepath, 'regopts.mat'), 'opts');
fpathsave       = fullfile(opts.savepath, 'reg_volume.tif');
if exist(fpathsave, 'file')
    delete(fpathsave);
end
options.compress = 'no';
saveastiff(newvol, fpathsave, options);
%--------------------------------------------------------------------------
end