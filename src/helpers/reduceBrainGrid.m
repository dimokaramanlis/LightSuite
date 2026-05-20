function bgdatared = reduceBrainGrid(bgdata, nred)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


segbins = [0;find(all(bgdata==0,2))];
finbins = cell(numel(segbins)-1,1);
for ibin = 1:numel(segbins)-1
    bgcurr = bgdata(segbins(ibin)+1:segbins(ibin+1)-1, :);
    Ncurr  = size(bgcurr, 1);
    samplevec =  linspace(1,Ncurr, round(Ncurr/nred));
    if Ncurr > 1
        bgcurr = interp1(1:Ncurr, single(bgcurr), samplevec);
    end
    finbins{ibin} = [bgcurr; 0 0 0];
end
bgdatared = cat(1, finbins{:});

end