function volpath = resolveTracingVolumePath(savepath, regopts, channel, volumepath)
% RESOLVETRACINGVOLUMEPATH  Choose which 20 µm registration volume to display
%   for probe / fiber annotation.
%
%   VOLPATH = resolveTracingVolumePath(SAVEPATH, REGOPTS, CHANNEL, VOLUMEPATH)
%
%   Implant tracing does not have to use the same channel as registration:
%   every channel is downsampled onto the *same* 20 µm registration grid, so
%   points clicked on any channel's volume map to the atlas through the same
%   transform.  This picks the channel volume the annotation GUIs load for
%   display, leaving the coordinate transform unchanged.
%
%   Precedence:
%     1. VOLUMEPATH – explicit file (absolute, or a name/relative path
%                     resolved inside SAVEPATH). Use for non-standard names.
%     2. CHANNEL    – channel index; builds
%                     chan_<CHANNEL>_sample_register_<registres>um.tif.
%     3. (default)  – REGOPTS.regvolpath, i.e. the registration channel.
%
%   Errors if the resolved file does not exist.

    if nargin < 3; channel    = []; end
    if nargin < 4; volumepath = ''; end

    if ~isempty(volumepath)
        volumepath = char(volumepath);
        if exist(volumepath, 'file')
            volpath = volumepath;                       % absolute / valid path
        else
            [~, n, e] = fileparts(volumepath);          % resolve inside savepath
            volpath = fullfile(savepath, [n e]);
        end
    elseif ~isempty(channel)
        volpath = fullfile(savepath, ...
            sprintf('chan_%d_sample_register_%dum.tif', channel, regopts.registres));
    else
        [~, n, e] = fileparts(regopts.regvolpath);      % registration channel
        volpath = fullfile(savepath, [n e]);
    end

    if ~exist(volpath, 'file')
        error('resolveTracingVolumePath:notFound', ...
            ['Tracing volume not found:\n  %s\n' ...
             'Use the ''Channel'' (index) or ''VolumePath'' option to select ' ...
             'the channel volume to trace on.'], volpath);
    end
end
