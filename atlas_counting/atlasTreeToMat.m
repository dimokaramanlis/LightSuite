function atlasmat = atlasTreeToMat(sttree)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% thanks gemini

cellArray = sttree.structure_id_path;

Nrows = size(cellArray, 1);
numericArrays = cell(Nrows, 1);

for i = 1:Nrows
    str = cellArray{i};

    % Remove leading and trailing slashes
    str = regexprep(str, '^/', '');
    str = regexprep(str, '/$', '');

    % Split the string by '/'
    parts = strsplit(str, '/');

    % Convert the string parts to numbers
    numParts = numel(parts);
    numericArray = zeros(1, numParts); % Preallocate for efficiency
    for j = 1:numParts
        numericArray(j) = str2double(parts{j});
    end

    numericArrays{i} = numericArray;
end


atlasmat = CelltoMatUE(numericArrays);

end