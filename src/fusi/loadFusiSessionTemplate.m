function tempout = loadFusiSessionTemplate(pathcurr, finshape)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[pathout, filename, ext] = fileparts(pathcurr);
temppath = fullfile(pathout, sprintf('%s_template%s',filename, ext));
if exist(temppath, 'file')
    tempout = load(temppath);
    tempout = tempout.tempout;
else
    % we build the template
    sessdata = load(pathcurr);
    Nframes  = size(sessdata.I, 3);
    isampinit = round(0.1*Nframes):round(0.9*Nframes);
    medall  = median(getRandomValsFromArray(sessdata.I,2e4));

    Xmat = reshape(sessdata.I,[],Nframes);
    %----------------------------------------------------------------------
    % chooose best template
    Npervol     = 50;
    Nrand       = 10;
    allcorrvals = nan(Nframes, Nrand);
    iframesrand = getRandomValsFromArray(isampinit, Npervol*Nrand);
    iframesrand = reshape(iframesrand, Npervol, Nrand);
    for ii = 1:Nrand
        iframescurr = iframesrand(:, ii);
        medvolcurr  = median(sessdata.I(:,:,iframescurr), 3)/medall;
        allcorrvals(:,ii) = corr(reshape(sessdata.I, [],Nframes), medvolcurr(:));
    end
    [~,imax] = max(median(allcorrvals,1));
    % tempuse  = median(sessdata.I(:,:,iframesrand(:, imax)), 3)/medall;
    % tempuse  = reshape(tempuse, finshape);
    %----------------------------------------------------------------------
    % refine template
    tempcorrs = allcorrvals(:, imax);
    ikeep     = find(tempcorrs > quantile(tempcorrs,0.9));
    tempout  = median(sessdata.I(:,:,ikeep), 3)/medall;
    tempout  = reshape(tempout, finshape);
    %----------------------------------------------------------------------
    save(temppath, 'tempout');
    %----------------------------------------------------------------------
end


end