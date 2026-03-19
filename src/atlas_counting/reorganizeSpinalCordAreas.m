function resout = reorganizeSpinalCordAreas(signals, volumes, parcelinfo, avinds, aggtype)
% REORGANIZESPINALCORDAREAS Aggregates spinal cord data to higher anatomical levels.
%
%   resout = reorganizeSpinalCordAreas(signals, volumes, parcelinfo, avinds, aggtype)
%
%   INPUTS:
%       signals    : Matrix of (Nareas x Nsegments x Nmice) containing intensity data.
%       volumes    : Matrix of (Nareas x Nsegments x Nmice) containing voxel counts.
%                    (If volumes is (Nareas x Nsegments), it broadcasts across mice).
%       parcelinfo : Table containing the Atlas_Regions.csv data.
%       avinds     : Vector of size (Nareas x 1) containing the Atlas IDs 
%                    corresponding to the rows of signals/volumes.
%       aggtype    : String, either 'division' (GM/WM) or 'structure' (Laminas+Funiculi).
%
%   OUTPUTS:
%       resout     : Struct containing .counts, .signal, .names, .indices
%
%   Dependencies: 
%       Requires parcelinfo to have columns: 'id', 'parent_ID', 'acronym', 'name'.

    %----------------------------------------------------------------------
    % 1. Handle Input Dimensions
    %----------------------------------------------------------------------
    [~, Nsegs, Nmice] = size(signals);
    
    % Handle case where volume is not per-mouse (common atlas volume)
    if ndims(volumes) == 2
        volumes = repmat(volumes, 1, 1, Nmice);
    end

    %----------------------------------------------------------------------
    % 2. Define Aggregation Targets based on Type
    %----------------------------------------------------------------------
    switch lower(aggtype)
        case 'division'
            % Targets: Gray Matter (GM), White Matter (WM)
            targetAcronyms = {'GM', 'WM'};
            
        case 'structure'
            % Targets: Laminas I-X Combined (201-210) + df, lf, vf
            % We search for specific IDs for Laminas and Acronyms for Funiculi
            
            % Laminas 1-10 Combined usually have IDs 201-210 in this specific atlas
            laminaIDs = 201:210; 
            % Funiculi Acronyms (Dorsal, Lateral, Ventral)
            funiculiAcronyms = {'df', 'lf', 'vf'};
            targetAcronyms = [arrayfun(@num2str, laminaIDs, 'UniformOutput', false), funiculiAcronyms];
        otherwise
            error('aggtype must be either "division" or "structure"');
    end

    %----------------------------------------------------------------------
    % 3. Initialize Outputs
    %----------------------------------------------------------------------
    Ntargets = numel(targetAcronyms);
    
    sigout   = nan(Ntargets, Nsegs, Nmice);
    volout   = nan(Ntargets, Nsegs, Nmice);
    nameout  = cell(Ntargets, 3); % Name, Acronym, Parent(placeholder)
    indsout  = nan(Ntargets, 1);

    %----------------------------------------------------------------------
    % 4. Main Aggregation Loop
    %----------------------------------------------------------------------
    for ii = 1:Ntargets
        currTarget = targetAcronyms{ii};
        
        % A. Find the ID of the target area in parcelinfo
        if all(isstrprop(currTarget, 'digit')) 
            % It's a numeric ID (like '201')
            targetID = str2double(currTarget);
            rowIdx   = find(parcelinfo.id == targetID, 1);
        else
            % It's an acronym (like 'GM' or 'df')
            rowIdx = find(strcmp(parcelinfo.acronym, currTarget), 1);
            if isempty(rowIdx) && strcmp(currTarget, 'lf')
                 % Fallback: sometimes 'lfc' (Lateral Funiculus Complete) is used if 'lf' missing
                 rowIdx = find(strcmp(parcelinfo.acronym, 'lfc'), 1);
            end
            if ~isempty(rowIdx)
                targetID = parcelinfo.id(rowIdx);
            else
                warning('Target %s not found in parcelinfo', currTarget);
                continue;
            end
        end
        
        if isempty(rowIdx)
            continue; 
        end
        
        % Store Metadata
        nameout(ii, :) = {parcelinfo.name{rowIdx}, parcelinfo.acronym{rowIdx}, aggtype};
        indsout(ii)    = targetID;

        % B. Find all raw Atlas IDs (avinds) that belong to this Target
        % This involves finding children, grandchildren, etc.
        relevantIDs = getDescendants(targetID, parcelinfo);
        
        % Include the target itself if it exists in the data
        relevantIDs = [relevantIDs; targetID]; 
        
        % Find which rows in our input matrices correspond to these IDs
        dataRows = find(ismember(avinds, relevantIDs));
        
        if isempty(dataRows)
            continue; % No data for this region
        end

        % C. Aggregate Volumes (Sum)
        % volumes(dataRows, :, :) -> Sum across rows
        vol_subset = volumes(dataRows, :, :);
        currVol    = sum(vol_subset, 1, 'omitmissing'); 
        volout(ii, :, :) = currVol;

        % D. Aggregate Signals (Weighted Average by Volume)
        % Weights = Volume of area / Total Volume of target
        % Note: We calculate weights per segment/mouse to handle missing data correctly
        
        sig_subset = signals(dataRows, :, :);
        
        % Avoid division by zero
        safeVol = currVol;
        safeVol(safeVol == 0) = 1; 
        
        % Weight calculation
        weights = vol_subset ./ safeVol;
        
        % Weighted sum
        weightedSig = sum(sig_subset .* weights, 1, 'omitmissing');
        sigout(ii, :, :) = weightedSig;
        
    end

    %----------------------------------------------------------------------
    % 5. Format Output
    %----------------------------------------------------------------------
    % Remove empty rows if any targets were totally missing from atlas definition
    keepIdx = ~cellfun('isempty', nameout(:,1));
    
    resout.counts   = volout(keepIdx, :, :); % "counts" in your original code referred to volume/cells
    resout.volumes  = volout(keepIdx, :, :);
    resout.signal   = sigout(keepIdx, :, :);
    resout.names    = nameout(keepIdx, :);
    resout.indices  = indsout(keepIdx);
end

%--------------------------------------------------------------------------
% Helper: Recursive Descendant Finder
%--------------------------------------------------------------------------
function allChildren = getDescendants(parentID, parcelinfo)
    % Finds all IDs that are children (direct or indirect) of parentID
    
    % 1. Find direct children based on parent_ID column
    directChildrenIdx = find(parcelinfo.parent_ID == parentID);
    directChildrenIDs = parcelinfo.id(directChildrenIdx);
    
    % 2. Special Case for "Combined" layers (200 series) in this specific atlas
    % These often list their constituent IDs in 'children_IDs' column if available,
    % or we infer them. If the CSV has 'children_IDs' as a string:
    if ismember('children_IDs', parcelinfo.Properties.VariableNames)
        pIdx = find(parcelinfo.id == parentID, 1);
        if ~isempty(pIdx)
            childStr = parcelinfo.children_IDs{pIdx};
            if ~isempty(childStr) && ischar(childStr)
                % Parse string "1, 2, 3" or "1Sp, 1"
                % Simple regex to grab digits
                parsedIDs = str2double(regexp(childStr, '\d+', 'match'));
                directChildrenIDs = unique([directChildrenIDs; parsedIDs(:)]);
            elseif isnumeric(childStr)
                 directChildrenIDs = unique([directChildrenIDs; childStr(:)]);
            end
        end
    end
    
    % 3. Recurse
    allChildren = directChildrenIDs;
    if ~isempty(directChildrenIDs)
        for k = 1:length(directChildrenIDs)
            grandChildren = getDescendants(directChildrenIDs(k), parcelinfo);
            allChildren = [allChildren; grandChildren]; %#ok<AGROW>
        end
    end
    allChildren = unique(allChildren);
end