function [tformout, res] = alignAtlasToSample(pcatlas, pcsamp, tformslices)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

rng(1);
Xlocs          = pcsamp.Location;
[apun, ~, ic]  = unique(pcsamp.Location(:, 2));
Nslices        = numel(tformslices);
Nptsdata       = size(Xlocs, 1);
Nptsatlas      = pcatlas.Count;
Nptsperslice   = ceil(10000/Nslices);
assert(numel(apun) == Nslices)

slicedist = median(diff(apun));

fprintf('Aligning atlas to 3D slice volume...\n'); tic;
% 
Xdown = cell(Nslices, 1);
for islice = 1:Nslices
    Xcurr = Xlocs(ic == islice, :);
    Ncurr  = size(Xcurr, 1);
    % careful with the dimensions here!!!!!!!
    Xcurr(:, [3 1]) = tformslices(islice).transformPointsForward(Xcurr(:, [3 1]));
    % Xcurr(:, 2)     = Xcurr(:, 2) + randn([Ncurr, 1]) * slicedist/4;
    % !!!
    pccurr = pointCloud(Xcurr);
    pccurr = pcdownsample(pccurr, 'random', Nptsperslice/Ncurr);

    Xcurr         = pccurr.Location;
    Xdown{islice} = Xcurr;

    % plot(pccurr.Location(:,3), pccurr.Location(:,1), '.',...
    %     pccurr2.Location(:,3), pccurr2.Location(:,1), '.')
    % ax = gca; ax.YDir = 'reverse';
    % pause;
end
Xdown   = cat(1, Xdown{:});
pcplot  = pointCloud(Xdown);

% Xdown = cell(Nslices, 1);
% for islice = 1:Nslices
%     Xcurr = Xlocs(ic == islice, :);
%     Xcurr(:, [1 3]) = tformslices(islice).transformPointsForward(Xcurr(:, [1 3]));
%     Xdown{islice} = Xcurr;
% end
% Xdown   = cat(1, Xdown{:});
% pcplot  = pointCloud(Xdown);
% pcplot =  pcdownsample(pcplot, 'nonuniformGridSample', ceil(Nptsdata/1e4), 'PreserveStructure',true);
% 

% pcplot = pcdownsample(pcsamp, 'random', 0.05, 'PreserveStructure',true);

% pcmov  = pcdownsample(pcatlas, 'random', 1e4/Nptsatlas, 'PreserveStructure',true);
pcmov  = pcdownsample(pcatlas, 'nonuniformGridSample',  ceil(Nptsatlas/0.5e4));

[tformout, pcreg, res] = pcregistercpd(pcmov, pcplot, "Transform","Rigid",...
    "Verbose",false,"OutlierRatio",0.00, 'MaxIterations', 10, 'Tolerance', 1e-6);

% figure;
% pcshowpair(pctransform(pcplot, tformout.invert), pcmov);

% pcshowpair(pcplot, pcmov);



fprintf('Done! Took %2.2f s. Overall 3d error = %2.3f\n', toc, res)

end