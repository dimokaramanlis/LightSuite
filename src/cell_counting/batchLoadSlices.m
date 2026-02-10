function slicevol = batchLoadSlices(slicepaths,sliceids, sliceyx)
%BATCHLOADSLICES Summary of this function goes here
%   Detailed explanation goes here

slicevol = zeros(sliceyx(1), sliceyx(2), numel(sliceids), 'uint16');
for ii = 1:numel(sliceids)
    islice = sliceids(ii);
    slicevol(:, :, ii) = imread(fullfile(slicepaths(islice).folder, slicepaths(islice).name));
end


end

