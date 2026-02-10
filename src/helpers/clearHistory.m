function clearHistory()
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    % remove old folders
olddirs = dir(fullfile(tempdir, 'transformix*'));
for idir = 1:numel(olddirs)
    rmdir(fullfile(olddirs(idir).folder, olddirs(idir).name), 's')
end

end