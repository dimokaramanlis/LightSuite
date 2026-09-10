function [areastat, areavols] = cordAreaStatistics(vol, av, grouping, areafun, atlasres)
%CORDAREASTATISTICS Summarise a registered cord volume per region and segment.
%
%   [AREASTAT, AREAVOLS] = CORDAREASTATISTICS(VOL, AV, GROUPING, AREAFUN,
%   ATLASRES) reduces the registered volume VOL to one number per atlas region
%   and per spinal segment, using AREAFUN (@median if not given). AV is the
%   annotation on the same grid as VOL and GROUPING comes from
%   CORDATLASGROUPING.
%
%   The cord equivalent of what GENERATEREGISTEREDBRAINVOLUMES does per
%   hemisphere: there is no left/right split here, because a cord is registered
%   as one object and the interesting axis is the rostrocaudal one instead.
%
%   AREASTAT is Nareas x Nsegments, NaN where a region does not reach a segment.
%   Region 0 - everything outside the cord - is kept as the background level,
%   which is what per-region intensities are normally expressed relative to.
%   AREAVOLS is the volume of every region in every segment, in mm^3.
%
%   The volume is walked one segment at a time so that no whole-volume index
%   array is ever built.
%
%   See also CORDATLASGROUPING, GENERATEREGISTEREDCORDVOLUME,
%   GENERATEREGISTEREDBRAINVOLUMES.

if nargin < 4 || isempty(areafun)
    areafun = @median;
end

Nareas    = grouping.Nareas;
Nsegments = grouping.Nsegments;

areastat = nan(Nareas, Nsegments, 'single');
areavols = zeros(Nareas, Nsegments, 'single');

voxvol = prod(atlasres) * 1e-9;   % um^3 -> mm^3

for iseg = 1:Nsegments
    iplanes = grouping.segofplane == iseg;
    if ~any(iplanes)
        continue
    end

    avblock  = av(:, :, iplanes);
    valblock = vol(:, :, iplanes);

    subs = double(grouping.lut(double(avblock(:)) + 1));
    % accumarray needs the fill value to match the class the function returns,
    % so everything goes through double and is narrowed afterwards
    vals = double(valblock(:));

    areastat(:, iseg) = single(accumarray(subs, vals, [Nareas 1], areafun, NaN));
    areavols(:, iseg) = single(accumarray(subs, 1,    [Nareas 1], @sum)) * voxvol;
end

end
