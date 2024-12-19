function cellData = loadColmTileInfo(xmlFile)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Initialize a structure array to store the cell data

xDoc = xmlread(xmlFile);
cells = xDoc.getElementsByTagName('Cell');
cellData = struct();

% Loop through each cell
for i = 0:cells.getLength-1
    currentCell = cells.item(i);
    
    % Extract ID
    id = str2double(currentCell.getElementsByTagName('ID').item(0).getFirstChild.getData);
    
    % Initialize structure for this cell
    cellData(i+1).ID = id;
    
    % Extract Position Data
    position = currentCell.getElementsByTagName('Position').item(0);
    x = str2double(position.getElementsByTagName('X').item(0).getFirstChild.getData);
    y = str2double(position.getElementsByTagName('Y').item(0).getFirstChild.getData);
    cellData(i+1).Position = [x, y];
    
    % Extract data for each HS_n (where n is variable)
    hsIndex = 0;
    while true
        hsTag = sprintf('HS_%d', hsIndex);
        hsNodeList = currentCell.getElementsByTagName(hsTag);
        if hsNodeList.getLength == 0
            break; % Exit if there are no more HS_n tags
        end
        hsNode = hsNodeList.item(0);
        
        % Extract parameters from each HS
        sample_mm = str2double(hsNode.getElementsByTagName('Sample_mm').item(0).getFirstChild.getData);
        detect_piezo_1_um = str2double(hsNode.getElementsByTagName('Detect_Piezo_1_um').item(0).getFirstChild.getData);
        detect_piezo_2_um = str2double(hsNode.getElementsByTagName('Detect_Piezo_2_um').item(0).getFirstChild.getData);
        hor_offset_1_um = str2double(hsNode.getElementsByTagName('Hor_Offset_1_um').item(0).getFirstChild.getData);
        hor_offset_2_um = str2double(hsNode.getElementsByTagName('Hor_Offset_2_um').item(0).getFirstChild.getData);
        laser_percent = str2double(hsNode.getElementsByTagName('Laser_Percent').item(0).getFirstChild.getData);
        
        % Store these details in the structure
        cellData(i+1).(hsTag) = struct('Sample_mm', sample_mm, ...
                                        'Detect_Piezo_1_um', detect_piezo_1_um, ...
                                        'Detect_Piezo_2_um', detect_piezo_2_um, ...
                                        'Hor_Offset_1_um', hor_offset_1_um, ...
                                        'Hor_Offset_2_um', hor_offset_2_um, ...
                                        'Laser_Percent', laser_percent);
        hsIndex = hsIndex + 1;
    end
end

end