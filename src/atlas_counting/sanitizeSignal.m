function sigout = sanitizeSignal(signalcheck, bkgthres)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[Nareas, ~, Nmice] = size(signalcheck);
sigout             = true(size(signalcheck));
Nsd = 2;
for imouse = 1:Nmice
    allsigs = signalcheck(:, :, imouse);
    meanval = median(allsigs, 'omitmissing');
    sdval   = robustStd(allsigs) * Nsd;
    inanbk    = (allsigs < (meanval - sdval)) | (allsigs < bkgthres);

    sigout(:, :, imouse) = ~inanbk;
    % inanplot = any(inan, 2);
    % alldiff = diff(allsigs, [],2);
    % maxc = max(allsigs, [],'all');
    % plot(allsigs(:,1),allsigs(:,2),'o',...
    %     allsigs(inanplot,1),allsigs(inanplot,2),'o', [0 maxc],[0 maxc]);
    % axis equal;
    % pause;
    % median(alldiff)
end
end