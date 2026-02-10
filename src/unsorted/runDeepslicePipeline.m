function [status, cmdout] = runDeepslicePipeline(pythonExecutable, pythonScriptPath, deepSliceParentDir, imageFolderPath, varargin)
%RUN_DEEPSLICE_PIPELINE Runs the DeepSlice Python processing script from MATLAB.
%
%   [status, cmdout] = RUN_DEEPSLICE_PIPELINE(pythonExecutable, ...
%                                           pythonScriptPath, ...
%                                           deepSliceParentDir, ...
%                                           imageFolderPath, ...
%                                           Name, Value, ...)
%
%   Inputs:
%       pythonExecutable     - Full path to the Python executable in the
%                              'deepslice' conda environment.
%                              (e.g., 'C:\Users\YourUser\miniconda3\envs\deepslice\python.exe')
%       pythonScriptPath     - Full path to the 'run_deepslice.py' script.
%                              (e.g., 'C:\path\to\scripts\run_deepslice.py')
%       deepSliceParentDir   - Full path to the main DeepSlice project directory
%                              (containing the 'DeepSlice' module).
%                              (e.g., 'D:\Tools\DeepSlice_Main_Folder')
%       imageFolderPath      - Full path to the folder containing images for DeepSlice.
%                              (e.g., 'J:\AM152\lightsuite\deepslice_input_images\')
%
%   Optional Name-Value Pair Arguments:
%       'Species'            - Species for the model. String.
%                              Default: 'mouse'. Options: 'mouse', 'rat'.
%       'UseEnsemble'        - Whether to use the ensemble model. Logical.
%                              Default: true.
%       'UseSectionNumbers'  - Whether section numbers are in filenames. Logical.
%                              Default: true.
%
%   Outputs:
%       status               - Exit status of the system command (0 for success).
%       cmdout               - Captured output (stdout/stderr) from the Python script.
%
%   Example Usage (uncomment and modify paths to use):
%   ------------------------username=getenv('USERNAME')----------------------------
%   pyExec = 'C:\Users\YourUser\miniconda3\envs\deepslice\python.exe'; % <<< UPDATE THIS
%   pyScript = 'C:\path\to\your\scripts\run_deepslice.py';           % <<< UPDATE THIS
%   dsParent = 'D:\Research\Tools\DeepSlice_Main_Folder';             % <<< UPDATE THIS
%   imgFolder = 'J:\AM152\lightsuite\deepslice_test_images\';         % <<< UPDATE THIS
%
%   [status, output] = run_deepslice_pipeline(pyExec, pyScript, dsParent, imgFolder, ...
%                                             'Species', 'mouse', ...
%                                             'UseEnsemble', true, ...
%                                             'UseSectionNumbers', true);
%   if status == 0
%       disp('DeepSlice pipeline completed successfully.');
%   else
%       disp('DeepSlice pipeline failed. See output for details:');
%       disp(output);
%   end
%   ----------------------------------------------------
pythonExecutable = fullfile('C:\Users\', username,'AppData', 'Local','miniconda3','envs','deepslice','python.exe');
% --- Input Parser for Optional Arguments ---
p = inputParser;
addRequired(p, 'pythonExecutable', @(x) ischar(x) || isstring(x));
addRequired(p, 'pythonScriptPath', @(x) ischar(x) || isstring(x));
addRequired(p, 'deepSliceParentDir', @(x) ischar(x) || isstring(x));
addRequired(p, 'imageFolderPath', @(x) ischar(x) || isstring(x));

addParameter(p, 'Species', 'mouse', @(x) ismember(lower(x), {'mouse', 'rat'}));
addParameter(p, 'UseEnsemble', true, @islogical);
addParameter(p, 'UseSectionNumbers', true, @islogical);

parse(p, pythonExecutable, pythonScriptPath, deepSliceParentDir, imageFolderPath, varargin{:});

% Assign parsed inputs
pythonExecutable = char(p.Results.pythonExecutable);
pythonScriptPath = char(p.Results.pythonScriptPath);
deepSliceParentDir = char(p.Results.deepSliceParentDir);
imageFolderPath = char(p.Results.imageFolderPath);
species = p.Results.Species;
useEnsemble = p.Results.UseEnsemble;
useSectionNumbers = p.Results.UseSectionNumbers;

% --- Script Pre-flight Checks ---
disp('--- MATLAB Function: run_deepslice_pipeline ---');
if isempty(pythonExecutable)
    error('MATLAB:DeepSlice:ConfigError', 'Python executable path (`pythonExecutable`) is not set.');
end
% Check if python executable exists (basic check)
if exist(pythonExecutable, 'file') ~= 2 && ~(ispc && exist([pythonExecutable '.exe'], 'file') == 2) && ~(~ispc && exist(pythonExecutable, 'file') == 7)
    warning('MATLAB:DeepSlice:PathWarning', 'Python executable might not be found at the specified path: %s. Please verify.', pythonExecutable);
end

if isempty(pythonScriptPath) || exist(pythonScriptPath, 'file') ~= 2
    error('MATLAB:DeepSlice:ConfigError', 'Python script `run_deepslice.py` not found at: %s. Ensure the path is correct and the file exists.', pythonScriptPath);
end

if isempty(deepSliceParentDir) || exist(deepSliceParentDir, 'dir') ~= 7
    error('MATLAB:DeepSlice:ConfigError', 'DeepSlice parent directory (`deepSliceParentDir`) is not set or not found: %s. Please verify the path.', deepSliceParentDir);
end

if isempty(imageFolderPath) || exist(imageFolderPath, 'dir') ~= 7
    error('MATLAB:DeepSlice:ConfigError', 'Image folder path (`imageFolderPath`) is not set or not found: %s. Please verify the path.', imageFolderPath);
end

% --- Configure MATLAB's Python Environment ---
disp('Configuring MATLAB Python environment...');
try
    currentPyEnv = pyenv;
    % Check if the executable path needs to be updated or if Python is not loaded
    if ~strcmpi(currentPyEnv.Executable, pythonExecutable) || ~strcmpi(currentPyEnv.Status, 'Loaded')
        disp(['Attempting to set Python environment to: ', pythonExecutable]);
        pyenv('Version', pythonExecutable);
        updatedPyEnv = pyenv; % Fetch updated environment status
        if strcmpi(updatedPyEnv.Status, 'Loaded')
            disp(['Python environment successfully set to: ', updatedPyEnv.Executable]);
        else
            warning('MATLAB:DeepSlice:PyEnvWarning', 'Failed to load Python environment. Status: %s. Check MATLAB''s Python configuration.', updatedPyEnv.Status);
            error('MATLAB:DeepSlice:PyEnvError', 'Could not set the Python environment. Please check the `pythonExecutable` path and your MATLAB Python setup.');
        end
    else
        disp(['MATLAB Python environment already correctly set to: ', currentPyEnv.Executable]);
    end
catch ME
    disp(ME.message)
    error('MATLAB:DeepSlice:PyEnvError', 'Failed to set Python environment using pyenv. Ensure MATLAB''s Python interface is configured and the `pythonExecutable` path is correct. Error: %s', ME.message);
end

% --- Construct the System Command to Run the Python Script ---
% Quote paths to handle spaces correctly.
commandParts = {};
commandParts{end+1} = sprintf('"%s"', pythonExecutable);
commandParts{end+1} = sprintf('"%s"', pythonScriptPath);
commandParts{end+1} = sprintf('"%s"', imageFolderPath);
commandParts{end+1} = sprintf('"%s"', deepSliceParentDir);

commandParts{end+1} = sprintf('--species %s', species);

if useEnsemble
    commandParts{end+1} = '--ensemble';
else
    commandParts{end+1} = '--no-ensemble';
end

if useSectionNumbers
    commandParts{end+1} = '--section_numbers';
else
    commandParts{end+1} = '--no-section_numbers';
end

fullCommand = strjoin(commandParts, ' ');

% --- Execute the Python Script ---
disp(' ');
disp('Running DeepSlice Python script. This may take some time...');
disp('Output from Python script will be displayed below:');
disp(['Executing command: ', fullCommand]);
disp('----------------------------------------------------');

% Use the system command to execute the Python script.
% The '-echo' flag (on Windows) or similar behavior on other OSes
% helps display the Python script's stdout/stderr in the MATLAB command window.
[status, cmdout] = system(fullCommand, '-echo');

disp('----------------------------------------------------');
% --- Check Execution Status ---
if status == 0
    disp('Python script executed successfully.');
else
    warning('MATLAB:DeepSlice:PythonExecutionFailed', 'Python script execution failed with status code: %d.', status);
    disp('Review the output above for error messages from the Python script.');
end
disp(' ');
disp('--- MATLAB Function Finished ---');

end % End of function run_deepslice_pipeline

% To make this file runnable as a script for testing (with user-defined paths):
% You would typically call the function from another script or the command window.
% Example call (ensure paths are correctly set before uncommenting and running):
%{
% --- BEGIN USER CONFIGURATION FOR TESTING ---
% 1. Path to the Python executable within your "deepslice" conda environment.
% pythonExecutable_test = 'C:\Users\YourUser\miniconda3\envs\deepslice\python.exe'; % <<< !!! UPDATE THIS PATH !!!
% pythonExecutable_test = '/Users/youruser/miniconda3/envs/deepslice/bin/python';

% 2. Path to the 'run_deepslice.py' script.
%    If this MATLAB function and 'run_deepslice.py' are in the same directory,
%    this should work. Otherwise, provide the full path.
% pythonScriptPath_test = fullfile(fileparts(mfilename('fullpath')), 'run_deepslice.py');

% 3. Path to the main DeepSlice project directory.
% deepSliceParentDir_test = 'D:\Research\Tools\DeepSlice_Main_Folder'; % <<< !!! UPDATE THIS PATH !!!

% 4. Input folder path containing the images to be processed by DeepSlice.
% imageFolderPath_test = 'J:\AM152\lightsuite\deepslice_test_images\'; % <<< !!! UPDATE THIS PATH !!!

% 5. Optional Parameters for DeepSlice:
% species_test = 'mouse';         % 'mouse' or 'rat'
% useEnsemble_test = true;        % true or false
% useSectionNumbers_test = true;  % true or false
% --- END USER CONFIGURATION FOR TESTING ---

% if exist('pythonExecutable_test', 'var') % Check if test variables are set
%     disp('Running test call to run_deepslice_pipeline...');
%     [test_status, test_cmdout] = run_deepslice_pipeline(...
%         pythonExecutable_test, ...
%         pythonScriptPath_test, ...
%         deepSliceParentDir_test, ...
%         imageFolderPath_test, ...
%         'Species', species_test, ...
%         'UseEnsemble', useEnsemble_test, ...
%         'UseSectionNumbers', useSectionNumbers_test);
% 
%     if test_status == 0
%         fprintf('Test call successful. Output:\n%s\n', test_cmdout);
%     else
%         fprintf('Test call failed. Status: %d. Output:\n%s\n', test_status, test_cmdout);
%     end
% else
%     disp('Test variables (e.g., pythonExecutable_test) are not defined.');
%     disp('To test, uncomment and configure the test block above, then run this .m file.');
% end
%}
