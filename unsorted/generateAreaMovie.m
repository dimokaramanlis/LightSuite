function generateAreaMovie(av, parcelinfo, namesuse, valsuse, rotationAxis, Nframes, outputFilename)
    % --- Inputs ---
    if nargin < 4 || isempty(rotationAxis), rotationAxis = 'Z'; end % Default axis
    if nargin < 5 || isempty(Nframes), Nframes = 100; end           % Default frames
    if nargin < 6 || isempty(outputFilename)                       % Default filename/type
        [~, ~, ext] = fileparts(outputFilename);
        if isempty(ext)
           outputFilename = 'brain_rotation.mp4'; % Default to MP4
        end
    end

    % --- Determine Output Format ---
    [~, ~, ext] = fileparts(outputFilename);
    isGIF = strcmpi(ext, '.gif');
    isMP4 = strcmpi(ext, '.mp4');
    if ~isGIF && ~isMP4
        warning('Output format not supported or recognized. Defaulting to MP4.');
        outputFilename = [outputFilename '.mp4'];
        isMP4 = true;
        isGIF = false;
    end

    %======================================================================
    % --- Data Preparation ---
    % lowval         = 100;
    % areaskeep      = contains(parcelinfo.parcellation_term_name,  namesuse(:, 1)); % Assuming namesuse is column
    % idsfind        = unique(parcelinfo.parcellation_index(areaskeep));
    % 
    % avshow_full         = zeros(size(av), 'uint8');
    % avback_full         = zeros(size(av), 'uint8');
    % avback_full(av>0)   = lowval;
    % 
    % % Efficient highlighting using ismember
    % if ~isempty(idsfind)
    %     highlight_mask_full = ismembc(av, uint16(idsfind));
    %     avshow_full(highlight_mask_full) = 1; % Use 1 for overlay mask
    % end

    lowval              = 100;
    avback_full         = zeros(size(av), 'uint8');
    avback_full(av>0)   = lowval;
    avshow_full         = zeros(size(av), 'uint8');

    valscale = uint8(1 + 254*(valsuse-min(valsuse))/range(valsuse));
    tic;
    for ii = 1:numel(valsuse)
        areaskeep   = contains(parcelinfo.parcellation_term_name,  namesuse(ii, 1)); % Assuming namesuse is column
        idsfind     = unique(parcelinfo.parcellation_index(areaskeep));
        highlight_mask_full = ismembc(av, uint16(idsfind));
        avshow_full(highlight_mask_full) = valscale(ii); % Use 1 for overlay mask
    end
    toc;
   

    % % Efficient highlighting using ismember
    % if ~isempty(idsfind)
    %     highlight_mask_full = ismembc(av, uint16(idsfind));
    %     avshow_full(highlight_mask_full) = 1; % Use 1 for overlay mask
    % end

    % Downsampling
    resizeFactor = 0.5;
    avshow = imresize3(avshow_full, resizeFactor, "Method", "nearest");
    avback = imresize3(avback_full, resizeFactor, "Method", "nearest");
    clear avshow_full avback_full highlight_mask_full; % Free memory
    %======================================================================
    % --- Viewer Setup (User's Code) ---
    figview = uifigure('Name', 'Brain Rotation Movie', 'Visible', 'on'); % Start invisible
    figview.Units = 'centimeters';
    figview.Position = [1 1 25 20];
    viewerLabels = viewer3d(figview, BackgroundColor="white", BackgroundGradient="off", CameraZoom=1);
    hVol = volshow(avback, Parent=viewerLabels, RenderingStyle="GradientOpacity", OverlayData=avshow,...
        OverlayDataLimits=[0 255]);

    % Apply transparency settings
    hVol.Alphamap = linspace(0, 0.2, 256); % For avback (GradientOpacity uses this for gradient)
    hVol.OverlayAlphamap = 0.7;           % For avshow overlay

    % Define overlay colormap: Intensity 1 in avshow should be red
    overlay_cmap = inferno; % Red color for value 1
    hVol.OverlayColormap = overlay_cmap; % Apply the colormap

    viewerLabels.Toolbar = 'off'; % Turn off toolbar as requested

    % --- Camera and Movie Setup ---
    sz = size(avback);
    center = sz / 2 + 0.5; % Center of the volume
    viewerLabels.CameraTarget = center;

    % Determine camera distance based on volume size
    diagonal = norm(sz);
    camDistance = diagonal * 1.7; % Adjust multiplier as needed

    % Define angles for rotation
    angles = linspace(-pi, pi, Nframes+1);
    angles = angles(1:Nframes); % Use Nframes points

    % --- Setup Video Writer or GIF ---
    if isMP4
        writerObj = VideoWriter(outputFilename, 'MPEG-4');
        writerObj.FrameRate = min(30, Nframes/10); % Aim for ~10 sec movie
        if writerObj.FrameRate < 5, writerObj.FrameRate = 5; end
        open(writerObj);
        fprintf('Starting MP4 generation: %s\n', outputFilename);
    elseif isGIF
        gifDelayTime = 1 / min(30, Nframes/6); % Match approx framerate
         if gifDelayTime > 0.2, gifDelayTime = 0.2; end % Max delay 0.2s (5fps)
         if gifDelayTime < 0.02, gifDelayTime = 0.02; end % Min delay 0.02s (50fps)
        fprintf('Starting GIF generation: %s\n', outputFilename);
    end
%%
    % Make figure visible now that setup is done
    figview.Visible = 'on';
    pause(1.0); % Allow figure to render completely before first frame capture

    % --- Rotation and Frame Capture Loop ---
    camUpVector = [-1 0 0]; % Y is typically up

    for idx = 1:Nframes
        current_angle = angles(idx);
        relativePos = zeros(1, 3); % Initialize relative position vector

        % Calculate Camera Position based on rotationAxis
        switch upper(rotationAxis)
            case 'Z' % Rotation in XY plane
                relativePos(1) = camDistance * cos(current_angle); % X
                relativePos(2) = camDistance * sin(current_angle); % Y
                relativePos(3) = 0; % Keep Z slightly offset from center Z or fixed (e.g. 0 relative)
            case 'Y' % Rotation in XZ plane
                relativePos(1) = camDistance * cos(current_angle); % X
                relativePos(2) = 0; % Keep Y slightly offset from center Y or fixed (e.g. 0 relative)
                relativePos(3) = camDistance * sin(current_angle); % Z
            case 'X' % Rotation in YZ plane
                relativePos(1) = 0; % Keep X slightly offset from center X or fixed (e.g. 0 relative)
                relativePos(2) = camDistance * cos(current_angle); % Y
                relativePos(3) = camDistance * sin(current_angle); % Z
            otherwise
                error('Invalid rotationAxis. Choose ''X'', ''Y'', or ''Z''.');
        end


        % Set camera properties for this frame
        viewerLabels.CameraPosition = center + relativePos;
        viewerLabels.CameraUpVector = camUpVector;
        viewerLabels.CameraTarget = center; % Ensure target remains fixed

        drawnow; % Update the figure window

        % --- Capture Frame ---
    
        frame =  getframe(figview);


        % --- Write Frame to File ---
        if idx==1
            % don nothing
        else
            if isMP4
                writeVideo(writerObj, frame);
            elseif isGIF
                [indI, cm] = rgb2ind(frame, 256);
                if idx == 2
                    imwrite(indI, cm, outputFilename, "gif", LoopCount=inf, DelayTime=gifDelayTime);
                else
                    imwrite(indI, cm, outputFilename, "gif", WriteMode="append", DelayTime=gifDelayTime);
                end
            end
        end

        % Progress update
        if mod(idx, round(Nframes/10)) == 0 || idx == Nframes
             fprintf('Processed frame %d of %d (%.0f%%)\n', idx, Nframes, (idx/Nframes)*100);
        end
    end
%%
    % --- Finalize Video/GIF ---
    if isMP4
        close(writerObj);
        fprintf('MP4 generation complete: %s\n', outputFilename);
    elseif isGIF
         fprintf('GIF generation complete: %s\n', outputFilename);
    end

    % --- Cleanup ---
    % Optional: Close the figure window
    % close(figview);

end % End of function
