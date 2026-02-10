function volumefin = applyConsecutiveTransforms(volregister, ...
    volumeindices, transformarray, fillvalues)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[Ny, Nx, Nchan, ~] = size(volregister);
assert(numel(fillvalues) == Nchan)
raref     = imref2d([Ny, Nx]);
Nalign    = numel(volumeindices);

volumefin = volregister(:, :, :, volumeindices);

Tmat   = transformarray(1).A;
for ii = 1:Nalign
    %----------------------------------------------------------------------
    idx       = volumeindices(ii);
    Tmat      = Tmat*transformarray(ii).A;
    tformcurr = affinetform2d(Tmat);
    %----------------------------------------------------------------------
    for ichan = 1:Nchan
        currimage = volregister(:, :, ichan, idx);
        warpedim = imwarp(currimage, tformcurr, 'linear', ...
            'OutputView',raref, 'FillValues', fillvalues(ichan));
        volumefin(:, :, ichan, ii) = warpedim;
    end
    %----------------------------------------------------------------------
end

end