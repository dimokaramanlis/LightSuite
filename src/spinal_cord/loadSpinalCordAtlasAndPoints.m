function [tv, av, tvpts, atlasres, segmentinfo] = loadSpinalCordAtlasAndPoints(outputres)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%-------------------------------------------------------------------------
dplook   = fileparts(which("Segments.csv"));
dpav     = fullfile(dplook, "Annotation.tif");
dptv     = fullfile(dplook, "Template.tif");
tv       = tiffreadVolume(dptv);
av       = tiffreadVolume(dpav);
atlasres = [10 10 20];
%-------------------------------------------------------------------------
parcelinfo  = readtable(fullfile(dplook,'Atlas_Regions.csv'));

% parcelinfo  = readtable(fullfile(dplook,'structures.csv'));
% iwm         = strcmpi(parcelinfo.name, 'white matter');
iwm         = contains(lower(parcelinfo.name), 'fasciculus') | ...
    contains(lower(parcelinfo.name), 'funiculus');

wmchildren  = parcelinfo.children_IDs(iwm);
wmchildren  = cellfun(@(x) split(x, ','), wmchildren, 'UniformOutput',false);
wmids       = str2double(cat(1, wmchildren{:}));
wmids(isnan(wmids)) = [];
%-------------------------------------------------------------------------
segmentinfo  = readtable(fullfile(dplook,'Segments.csv'));
%-------------------------------------------------------------------------
tv       = imresize3(tv, 'Scale', atlasres./outputres);
av       = imresize3(av, 'nearest', 'Scale', atlasres./outputres);
%-------------------------------------------------------------------------
Nslices  = size(tv, 3);
ptlist   = cell(Nslices, 1);
thresuse = 1;
for islice = 1:Nslices
    avcurr  = av(:,:,islice);

    %======================================================================
    % % find white matter transitions
    iswmtrans = imgradient(ismember(avcurr, wmids)) > 0;
    % find extenal points
    isextpt =imgradient(avcurr>0) > 0;
    % white matter points
    ipts       = find(isextpt | iswmtrans);
    [rr,cc]    = ind2sub(size(avcurr), ipts);
    ptslice    = [cc,rr,islice*ones(size(rr))];
    ptcurr    = pointCloud(ptslice);
    ptcurr    = pcdownsample(ptcurr, 'nonuniformGridSample', 6);
    %======================================================================
    % scurr = single(medfilt2(tv(:,:,islice), [5 5]));
    % tthres = quantile(scurr(avcurr>0),0.5);
    % % gout  = gcurr./imgaussfilt(scurr, 5,"Padding","symmetric");
    % ipts       = find(scurr<tthres & avcurr>0);
    % [rr,cc]    = ind2sub(size(scurr), ipts);
    % ptslice    = [cc,rr,islice*ones(size(rr))];
    % ptcurr    = pointCloud(ptslice);
    % ptcurr    = pcdownsample(ptcurr, 'nonuniformGridSample', 12);
    % 
    % ptlist{islice} = ptcurr.Location;   
    %======================================================================
    % scurr = single(medfilt2(tv(:,:,islice), [9 9]));
    % % scurr = single(imgaussfilt(tv(:,:,islice), 2));
    % gcurr = imgradient(scurr);
    % % gout  = gcurr./imgaussfilt(scurr, 5,"Padding","symmetric");
    % gout  = gcurr./scurr;
    % 
    % ipts       = find(gout>thresuse);
    % [rr,cc]    = ind2sub(size(scurr), ipts);
    % ptslice    = [cc,rr,islice*ones(size(rr))];
    % ptcurr    = pointCloud(ptslice);
    % ptcurr    = pcdownsample(ptcurr, 'nonuniformGridSample', 6);

    %======================================================================
    ptlist{islice} = ptcurr.Location;   
    % plot(ptcurr.Location(:,1), ptcurr.Location(:,2),'o')
    %======================================================================
end
%%
tvpts    = single(cat(1, ptlist{:}));
%-------------------------------------------------------------------------
end