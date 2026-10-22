function report = updateLightSuite(varargin)
%UPDATELIGHTSUITE Pull the latest changes for LightSuite and its git dependencies.
%
%   updateLightSuite() finds LightSuite and the dependency repositories on
%   the MATLAB path (matlab_elastix, yamlmatlab, bcpd, Bioformats Image
%   Toolbox), fetches each one, prints the incoming commits and
%   fast-forwards it to its upstream branch.
%
%   updateLightSuite('checkonly', true) only fetches and lists the incoming
%   commits, without changing any files.
%
%   report = updateLightSuite(...) also returns a struct array with one
%   entry per repository (name, path, branch, status, commits).
%
%   Safety: only fast-forward pulls are done. A repository with diverged
%   history, or local edits that would be overwritten, is left untouched
%   and reported, so local work is never lost.
%
%   updateLightSuite('gitpath', 'C:\path\to\git.exe') uses that git.
%
%   Requires git. If it is not on the system PATH, the usual install
%   locations are searched: Git for Windows / Git Bash (registry and
%   Program Files), GitHub Desktop's bundled git, SourceTree, Scoop,
%   Chocolatey, Visual Studio, and on macOS/Linux Homebrew, Xcode tools and
%   GitHub Desktop.

p = inputParser;
p.addParameter('checkonly', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('gitpath', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
checkonly = logical(p.Results.checkonly);

[gitcmd, gitsource] = findGit(char(p.Results.gitpath));
if isempty(gitcmd)
    error('updateLightSuite:nogit', ['git was not found on the system PATH or ' ...
        'in the usual install locations (Git for Windows, GitHub Desktop, ...).\n' ...
        'Install it from https://git-scm.com, or pass its location with ' ...
        'updateLightSuite(''gitpath'', ''C:\\path\\to\\git.exe'').']);
end
if ~strcmp(gitsource, 'system PATH')
    fprintf('Using git from %s: %s\n', gitsource, gitcmd);
end

% Never let git wait for a password prompt that MATLAB cannot answer
oldprompt = getenv('GIT_TERMINAL_PROMPT');
setenv('GIT_TERMINAL_PROMPT', '0');
restoreEnv = onCleanup(@() setenv('GIT_TERMINAL_PROMPT', oldprompt));

%--------------------------------------------------------------------------
% Repositories: name and a file that lives inside it
reponames   = {'LightSuite', 'matlab_elastix', 'yamlmatlab', 'bcpd', ...
               'Bioformats Image Toolbox'};
markerfiles = {[mfilename('fullpath') '.m'], 'invertElastixTransformCP.m', ...
               'yaml.ReadYaml', 'bcpd.exe', 'BioformatsImage.m'};

report = struct('name', {}, 'path', {}, 'branch', {}, 'status', {}, 'commits', {});
seenroots = {};

fprintf('\n%s\n', repmat('==', [1 30]));
if checkonly
    fprintf('Checking for LightSuite updates (no files will change)...\n');
else
    fprintf('Updating LightSuite and its dependencies...\n');
end
fprintf('%s\n', repmat('==', [1 30]));

for irepo = 1:numel(reponames)
    entry = struct('name', reponames{irepo}, 'path', '', 'branch', '', ...
                   'status', '', 'commits', {{}});

    % Locate the repository root
    if isfile(markerfiles{irepo})
        markerpath = markerfiles{irepo};
    else
        markerpath = which(markerfiles{irepo});
    end
    if isempty(markerpath)
        entry.status = 'not found on MATLAB path';
        report = printAndStore(report, entry);
        continue
    end
    [status, root] = runGit(gitcmd, fileparts(markerpath), 'rev-parse --show-toplevel');
    if status ~= 0
        entry.status = 'not a git repository (skipped)';
        entry.path   = fileparts(markerpath);
        report = printAndStore(report, entry);
        continue
    end
    entry.path = root;
    if any(strcmpi(seenroots, root))
        continue
    end
    seenroots{end+1} = root; %#ok<AGROW>

    [~, entry.branch] = runGit(gitcmd, root, 'rev-parse --abbrev-ref HEAD');
    if strcmp(entry.branch, 'HEAD')
        entry.status = 'detached HEAD (skipped)';
        report = printAndStore(report, entry);
        continue
    end
    [status, ~] = runGit(gitcmd, root, 'rev-parse --abbrev-ref --symbolic-full-name @{u}');
    if status ~= 0
        entry.status = 'no upstream branch (skipped)';
        report = printAndStore(report, entry);
        continue
    end

    % Fetch and list incoming commits
    [status, msg] = runGit(gitcmd, root, 'fetch --quiet');
    if status ~= 0
        entry.status = sprintf('fetch failed: %s', firstLine(msg));
        report = printAndStore(report, entry);
        continue
    end
    [~, nbehind] = runGit(gitcmd, root, 'rev-list --count --no-merges HEAD..@{u}');
    [~, nahead]  = runGit(gitcmd, root, 'rev-list --count @{u}..HEAD');
    nbehind = str2double(nbehind);
    nahead  = str2double(nahead);

    if nbehind == 0
        entry.status = 'up to date';
        report = printAndStore(report, entry);
        continue
    end

    [~, logtxt] = runGit(gitcmd, root, ...
        'log --no-merges --date=short --pretty=format:"%h  %ad  %s (%an)" HEAD..@{u}');
    entry.commits = splitlines(string(logtxt));
    entry.commits = cellstr(entry.commits(strlength(entry.commits) > 0));

    if checkonly
        entry.status = sprintf('%d update(s) available', nbehind);
    elseif nahead > 0
        entry.status = sprintf(['%d update(s) available, but you have %d local ' ...
            'commit(s) - pull manually'], nbehind, nahead);
    else
        [status, msg] = runGit(gitcmd, root, 'pull --ff-only --quiet');
        if status == 0
            entry.status = sprintf('updated with %d new commit(s)', nbehind);
        else
            entry.status = sprintf(['%d update(s) available, but pull failed ' ...
                '(local edits in the way?): %s'], nbehind, firstLine(msg));
        end
    end
    report = printAndStore(report, entry);
end

%--------------------------------------------------------------------------
% Summary
nupdated = sum(startsWith({report.status}, 'updated'));
navail   = sum(contains({report.status}, 'available'));
fprintf('%s\n', repmat('==', [1 30]));
if nupdated > 0
    rehash path;
    fprintf('%d repository(ies) updated. Restart running GUIs to use the new code.\n', nupdated);
end
if navail > 0
    fprintf('%d repository(ies) have updates that were not pulled.\n', navail);
end
if nupdated == 0 && navail == 0
    fprintf('Everything is up to date.\n');
end
fprintf('%s\n\n', repmat('==', [1 30]));

if nargout == 0
    clear report
end
end

%==========================================================================
function [status, out] = runGit(gitcmd, repodir, args)
[status, out] = system(sprintf('%s -C "%s" %s', gitcmd, repodir, args));
out = strtrim(out);
end

function [gitcmd, source] = findGit(userpath)
% Return a command that runs git (quoted full path, or 'git' if on PATH),
% plus a short description of where it was found. Empty if none works.
gitcmd = '';
source = '';

% 1. Location given by the user
if ~isempty(userpath)
    if isfolder(userpath)
        userpath = fullfile(userpath, gitExeName());
    end
    if isWorkingGit(userpath)
        gitcmd = quotePath(userpath); source = 'user-supplied path';
    end
    return
end

% 2. System PATH
if isWorkingGit('git')
    gitcmd = 'git'; source = 'system PATH';
    return
end

% 3. Known install locations
candidates = {};   % {path, description}
if ispc
    progfiles = {getenv('ProgramFiles'), getenv('ProgramW6432'), ...
                 getenv('ProgramFiles(x86)'), 'C:\Program Files'};
    progfiles = unique(progfiles(~cellfun(@isempty, progfiles)), 'stable');
    localapp  = getenv('LOCALAPPDATA');
    userprof  = getenv('USERPROFILE');

    % Git for Windows / Git Bash: install folder from the registry first
    for hive = {'HKEY_LOCAL_MACHINE', 'HKEY_CURRENT_USER'}
        try %#ok<TRYNC>
            installdir = winqueryreg(hive{1}, 'SOFTWARE\GitForWindows', 'InstallPath');
            candidates(end+1, :) = {fullfile(installdir, 'cmd', 'git.exe'), 'Git for Windows'}; %#ok<AGROW>
        end
    end
    for ip = 1:numel(progfiles)
        candidates(end+1, :) = {fullfile(progfiles{ip}, 'Git', 'cmd', 'git.exe'), 'Git for Windows'}; %#ok<AGROW>
        candidates(end+1, :) = {fullfile(progfiles{ip}, 'Git', 'bin', 'git.exe'), 'Git Bash'}; %#ok<AGROW>
    end
    if ~isempty(localapp)
        candidates(end+1, :) = {fullfile(localapp, 'Programs', 'Git', 'cmd', 'git.exe'), 'Git for Windows (user install)'};
        % GitHub Desktop ships its own git in app-<version> folders: newest first
        appdirs = dir(fullfile(localapp, 'GitHubDesktop', 'app-*'));
        [~, order] = sort([appdirs.datenum], 'descend');
        for ia = order
            candidates(end+1, :) = {fullfile(appdirs(ia).folder, appdirs(ia).name, ...
                'resources', 'app', 'git', 'cmd', 'git.exe'), 'GitHub Desktop'}; %#ok<AGROW>
        end
        candidates(end+1, :) = {fullfile(localapp, 'Atlassian', 'SourceTree', 'git_local', 'bin', 'git.exe'), 'SourceTree'};
    end
    if ~isempty(userprof)
        candidates(end+1, :) = {fullfile(userprof, 'scoop', 'shims', 'git.exe'), 'Scoop'};
    end
    candidates(end+1, :) = {fullfile(getenv('ProgramData'), 'chocolatey', 'bin', 'git.exe'), 'Chocolatey'};
    % Git bundled with Visual Studio
    for ip = 1:numel(progfiles)
        vsgits = dir(fullfile(progfiles{ip}, 'Microsoft Visual Studio', '*', '*', 'Common7', ...
            'IDE', 'CommonExtensions', 'Microsoft', 'TeamFoundation', 'Team Explorer', ...
            'Git', 'cmd', 'git.exe'));
        for iv = 1:numel(vsgits)
            candidates(end+1, :) = {fullfile(vsgits(iv).folder, vsgits(iv).name), 'Visual Studio'}; %#ok<AGROW>
        end
    end
else
    % MATLAB started from the Dock/launcher often misses the shell PATH
    candidates = {'/usr/bin/git', 'system'; '/usr/local/bin/git', 'system'; ...
                  '/opt/homebrew/bin/git', 'Homebrew'; ...
                  '/Library/Developer/CommandLineTools/usr/bin/git', 'Xcode tools'; ...
                  '/Applications/GitHub Desktop.app/Contents/Resources/app/git/bin/git', 'GitHub Desktop'};
end

for ic = 1:size(candidates, 1)
    if isfile(candidates{ic, 1}) && isWorkingGit(candidates{ic, 1})
        gitcmd = quotePath(candidates{ic, 1});
        source = candidates{ic, 2};
        return
    end
end
end

function ok = isWorkingGit(exe)
if ~strcmp(exe, 'git')
    exe = quotePath(exe);
end
[status, out] = system(sprintf('%s --version', exe));
ok = status == 0 && contains(out, 'git version');
end

function q = quotePath(p)
q = ['"' p '"'];
end

function name = gitExeName()
if ispc
    name = 'git.exe';
else
    name = 'git';
end
end

function s = firstLine(txt)
lines = splitlines(string(txt));
lines = lines(strlength(strtrim(lines)) > 0);
if isempty(lines)
    s = '';
else
    s = char(lines(1));
end
end

function report = printAndStore(report, entry)
if isempty(entry.branch)
    fprintf('%-26s %s\n', entry.name, entry.status);
else
    fprintf('%-26s [%s] %s\n', entry.name, entry.branch, entry.status);
end
for ic = 1:numel(entry.commits)
    fprintf('      %s\n', entry.commits{ic});
end
report(end+1) = entry;
end
