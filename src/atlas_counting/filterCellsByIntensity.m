function igood = filterCellsByIntensity(cellcoords, annvol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

finalpts = [round(cellcoords(:, 1:3)) cellcoords(:, 4:end)]; % rounding to transfer correspondance to pixels
%--------------------------------------------------------------------------
irem0 = any(finalpts(:,1:3)<1, 2);
iremx = finalpts(:,1) > size(annvol, 2);
iremy = finalpts(:,2) > size(annvol, 1);
iremz = finalpts(:,3) > size(annvol, 3);
irem  = irem0 | iremx | iremy | iremz;
%--------------------------------------------------------------------------
% intensity filtering
indcells     = sub2ind(size(annvol), finalpts(~irem,2), finalpts(~irem,1), finalpts(~irem,3));
irem(~irem)  = annvol(indcells) == 0;

mingoodT = quantile(finalpts(~irem,4), 0.01);
maxbadT  = mode(finalpts(irem,4));
if maxbadT > mingoodT
    maxbadT = 1;
end
Tfin     = (maxbadT + mingoodT)/2;
Tfin     = max(Tfin, 1);
    
% clf; hold on;
% histogram(finalpts(irem,4), linspace(0, 2e4), 'Normalization','probability')
% histogram(finalpts(~irem,4),linspace(0, 2e4), 'Normalization','probability')
% xline(Tfin, 'LineWidth',2)
%%
irem  = irem | (finalpts(:, 4) < Tfin);
igood = ~irem;
end