%==========================================================================
% 1. Check MATLAB requirements
fprintf('\n%s\n', repmat('==', [1 25]));
releaseInfo  = ver;
tf           = isMATLABReleaseOlderThan("R2022b");
if tf
    error("Your MATLAB version is older that R2022b");
else
    fprintf('Your MATLAB version is valid\n');
end
fprintf('Looking for required MATLAB toolboxes... ');
toolboxNames = {"Computer Vision Toolbox", ...
                "Image Processing Toolbox",...
                "Optimization Toolbox", ...
                "Parallel Computing Toolbox", ...
                "Statistics and Machine Learning Toolbox"};
for itoolbox = 1:numel(toolboxNames)
    hastoolbox = any(strcmp(toolboxNames{itoolbox}, {releaseInfo.Name}));
    if ~hastoolbox
        error("You're missing the %s", toolboxNames{itoolbox});
    end
end
fprintf('Done!\n');
%==========================================================================
% 2. Check repository requirements
fprintf('%s\n', repmat('==', [1 25]));
reponames        = {'matlab_elastix',                 'bcpd'};
reporelatedfiles = {'invertElastixTransformCP.m', 'bcpd.exe'};
for irepofile = 1:numel(reporelatedfiles)
    hasrepofile = which(reporelatedfiles{irepofile});
    if ~isempty(hasrepofile)
        melpath = fileparts(hasrepofile);
        fprintf('Found %s in %s\n', reponames{irepofile}, melpath);
    else
        error("You have to add %s available to MATLAB's path", reponames{irepofile})
    end
end
%==========================================================================
% 3. Check elastix requirements
fprintf('%s\n', repmat('==', [1 25]));
fprintf('Looking for the elastix/transformix binaries...\n');
commnames = {'elastix', 'transformix'};
for ii = 1:numel(commsuse)
    [status, result]=system(sprintf("%s --help", commnames{ii}));
    if status == 0

        strgaps = strfind(result, newline);
        disp(result(1:strgaps(1)-1))
    else
        error("You have to add %s to your system's path.", commnames{ii})
    end
end

fprintf('Found! Make sure the version is 5.1.0!\n');
%==========================================================================
% 4. Check atlases


fprintf('%s\n', repmat('==', [1 25]));
reponames        = {'Allen brain',                 'spinal cord'};
reporelatedfiles = {'average_template_10.nii.gz', 'Segments.csv'};
for irepofile = 1:numel(reporelatedfiles)
    hasrepofile = which(reporelatedfiles{irepofile});
    if ~isempty(hasrepofile)
        melpath = fileparts(hasrepofile);
        fprintf('Found %s atlas in %s\n', reponames{irepofile}, melpath);
    else
        error("The %s atlas was not found, maybe you'll run into problems", reponames{irepofile})
    end
end

%==========================================================================