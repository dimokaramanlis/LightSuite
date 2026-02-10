function [final_pts, tform_array] = iterative_straightening(samppts, iax, Nslices)
    % MAIN SCRIPT TO RUN THE ITERATIONS
    % samppts: [N x 3] point cloud
    % iax:     Index of the slice axis (e.g., 3 for Z)
    
    % Configuration
    N_iters = 5;
    batch_size = 50; % Adjust based on data density
    
    current_pts = samppts;
    
    % Initialize "Global" transform per slice as Identity
    % We track this to give you the final cumulative result
    final_tforms(Nslices, 1) = rigidtform2d; 
    
    fprintf('=== STARTING ROBUST STRAIGHTENING ===\n');
    
    for iter = 1:N_iters
        fprintf('\n--- Iteration %d / %d ---\n', iter, N_iters);
        
        % 1. MEASURE STRAIGHTNESS (Metric)
        % We use the smooth vector solver just to MEASURE how wavy it is.
        % Ideally, the offsets should go to zero.
        center_offsets = measure_straightness(current_pts, iax, Nslices);
        metric = norm(center_offsets(:)) / sqrt(Nslices);
        fprintf('   Current Straightness Error (RMS): %.4f\n', metric);
        
        % 2. COMPUTE UPDATE (Middle-Out Registration)
        % This returns the incremental transforms for this step
        [step_tforms, current_pts] = alignConsecutiveCordSlices(current_pts, iax, batch_size);
        
        % 3. ACCUMULATE TRANSFORMS
        % T_final = T_step * T_old
        for i = 1:Nslices
            final_tforms(i) = rigidtform2d(step_tforms(i).A * final_tforms(i).A);
        end
        
        % Break if converged
        if metric < 0.05 % Threshold depends on your units
            fprintf('   Converged!\n');
            break;
        end
    end
    
    final_pts = current_pts;
    tform_array = final_tforms;
    
    % Plot result
    figure; 
    subplot(1,2,1); pcshow(samppts); title('Original'); axis on;
    subplot(1,2,2); pcshow(final_pts); title('Straightened'); axis on;
end

function [tform_full, aligned_pts] = alignConsecutiveCordSlices(pcsamp, iax, batchsize)
    % ROBUST MIDDLE-OUT REGISTRATION WITH INTERPOLATION
    bcpd = 'C:\Users\karamanl\Documents\GitHub\bcpd\win\bcpd.exe';
    
    Nslices = max(pcsamp(:, iax));
    other_dims = setdiff(1:3, iax);
    
    % 1. DEFINE BATCHES
    Nbatches = ceil(Nslices/batchsize);
    mid_batch = ceil(Nbatches/2);
    
    % We will store transforms at "Anchor Slices"
    anchors = zeros(Nbatches, 1);
    batch_tforms(Nbatches, 1) = rigidtform2d; 
    valid_batch_mask = false(Nbatches, 1); % Track which batches actually worked
    
    % Calculate anchor positions (approximate center of each batch)
    for b = 1:Nbatches
        istart = (b-1)*batchsize + 1;
        iend   = min(b*batchsize, Nslices);
        anchors(b) = floor((istart + iend)/2);
    end
    
    aligned_pts = pcsamp;
    
    fprintf('   Aligning %d batches (Middle: %d)...\n', Nbatches, mid_batch);
    
    % =====================================================================
    % 2. FORWARD PASS (Middle -> End)
    % =====================================================================
    valid_batch_mask(mid_batch) = true; % Middle is Identity (Reference)
    
    for b = mid_batch+1 : Nbatches
        % Fixed: Batch b-1 (Previous neighbor)
        % Moving: Batch b (Current)
        [pcmov, pcfix] = get_batch_pts(aligned_pts, b, b-1, batchsize, Nslices, iax, other_dims);
        
        if isempty(pcmov) || isempty(pcfix)
            % If empty, we just copy the transform from the previous neighbor
            batch_tforms(b) = batch_tforms(b-1);
            valid_batch_mask(b) = true; 
            continue;
        end
        
        tform = run_bcpd(pcmov, pcfix, bcpd);
        
        % Update Points in place
        idx_start = (b-1)*batchsize + 1;
        idx_end   = min(b*batchsize, Nslices);
        mask = aligned_pts(:, iax) >= idx_start & aligned_pts(:, iax) <= idx_end;
        
        pts = aligned_pts(mask, other_dims);
        [u,v] = transformPointsForward(tform, pts(:,1), pts(:,2));
        aligned_pts(mask, other_dims) = [u,v];
        
        % Accumulate: T_current = T_rel * T_neighbor
        batch_tforms(b) = rigidtform2d(tform.A * batch_tforms(b-1).A);
        valid_batch_mask(b) = true;
    end
    
    % =====================================================================
    % 3. BACKWARD PASS (Middle -> Start)
    % =====================================================================
    for b = mid_batch-1 : -1 : 1
        % Fixed: Batch b+1 (Next neighbor)
        % Moving: Batch b (Current)
        [pcmov, pcfix] = get_batch_pts(aligned_pts, b, b+1, batchsize, Nslices, iax, other_dims);
        
        if isempty(pcmov) || isempty(pcfix)
             % Copy neighbor's transform if registration impossible
            batch_tforms(b) = batch_tforms(b+1);
            valid_batch_mask(b) = true;
            continue;
        end
        
        tform = run_bcpd(pcmov, pcfix, bcpd);
        
        idx_start = (b-1)*batchsize + 1;
        idx_end   = min(b*batchsize, Nslices);
        mask = aligned_pts(:, iax) >= idx_start & aligned_pts(:, iax) <= idx_end;
        
        pts = aligned_pts(mask, other_dims);
        [u,v] = transformPointsForward(tform, pts(:,1), pts(:,2));
        aligned_pts(mask, other_dims) = [u,v];
        
        % Accumulate: T_current = T_rel * T_neighbor
        batch_tforms(b) = rigidtform2d(tform.A * batch_tforms(b+1).A);
        valid_batch_mask(b) = true;
    end
    
    % =====================================================================
    % 4. ROBUST INTERPOLATION
    % =====================================================================
    fprintf('   Interpolating transforms for %d slices...\n', Nslices);
    tform_full(Nslices, 1) = rigidtform2d;
    
    % Only use valid anchors for interpolation
    valid_anchors = anchors(valid_batch_mask);
    valid_tforms  = batch_tforms(valid_batch_mask);
    
    % Safety check: If for some reason we have < 2 anchors, we can't interpolate.
    if numel(valid_anchors) < 2
        warning('Not enough valid batches for interpolation. Returning Identity.');
        return; 
    end

    % Extract arrays for interp1
    trans_vecs = zeros(numel(valid_anchors), 2);
    rot_angles = zeros(numel(valid_anchors), 1);
    
    for k = 1:numel(valid_anchors)
        T = valid_tforms(k).T; % [R 0; tx ty 1]
        trans_vecs(k, :) = T(3, 1:2);
        % Extract angle safely
        rot_angles(k) = atan2(T(1,2), T(1,1));
    end
    
    % Unwrap is critical for angles
    rot_angles = unwrap(rot_angles);
    
    % Clean duplicate anchors (rare, but causes interp1 crash)
    [valid_anchors, unique_idx] = unique(valid_anchors);
    trans_vecs = trans_vecs(unique_idx, :);
    rot_angles = rot_angles(unique_idx);

    % INTERPOLATE
    % We use 'linear' and 'extrap' to cover slices 1..Nslices
    xq = (1:Nslices)';
    
    theta_interp = interp1(valid_anchors, rot_angles, xq, 'linear', 'extrap');
    tx_interp    = interp1(valid_anchors, trans_vecs(:,1), xq, 'linear', 'extrap');
    ty_interp    = interp1(valid_anchors, trans_vecs(:,2), xq, 'linear', 'extrap');
    
    % Construct Transforms
    for i = 1:Nslices
        % 1. Get components
        th = theta_interp(i);
        tx = tx_interp(i);
        ty = ty_interp(i);
        
        % 2. Check for NaN (Stop the crash!)
        if isnan(th) || isnan(tx) || isnan(ty)
            th = 0; tx = 0; ty = 0; % Fallback to Identity
        end
        
        % 3. Explicit Rotation Matrix Construction
        % (Avoids rounding errors creating non-orthogonal matrices)
        c = cos(th); 
        s = sin(th);
        
        R = [c s; -s c]; % Orthogonal by definition
        
        % 4. Build 3x3
        T_mat = eye(3);
        T_mat(1:2, 1:2) = R;
        T_mat(3, 1:2) = [tx, ty];
        
        tform_full(i) = rigidtform2d(T_mat);
    end
    
    % Final clean application to ensure output points match the transforms
    aligned_pts = pcsamp; 
    for i = 1:Nslices
        mask = aligned_pts(:, iax) == i;
        if any(mask)
            pts = aligned_pts(mask, other_dims);
            [u, v] = transformPointsForward(tform_full(i), pts(:,1), pts(:,2));
            aligned_pts(mask, other_dims) = [u, v];
        end
    end
end

% --- HELPER: DATA EXTRACTION ---
function [pcmov, pcfix] = get_batch_pts(pts, b_mov, b_fix, bsize, N, iax, odims)
    s_mov = (b_mov-1)*bsize + 1; e_mov = min(b_mov*bsize, N);
    s_fix = (b_fix-1)*bsize + 1; e_fix = min(b_fix*bsize, N);
    
    mask_mov = pts(:, iax) >= s_mov & pts(:, iax) <= e_mov;
    mask_fix = pts(:, iax) >= s_fix & pts(:, iax) <= e_fix;
    
    pcmov = pts(mask_mov, odims);
    pcfix = pts(mask_fix, odims);
end

% --- HELPER: RUN BCPD ---
function tform = run_bcpd(mov, fix, bcpd_path)
    [~, bfit] = pcregisterBCPD(mov, fix, 'TransformType','Rigid',...
        'BCPDPath', bcpd_path, 'OutlierRatio', 0.1, ...
        'Gamma', 1, 'Verbose', false, 'ConvergenceTolerance', 1e-6);
    tform = affinetform2d(bfit);
end

% --- HELPER: MEASURE ERROR ---
function offsets = measure_straightness(pts, iax, N)
    % Simple logic: Calculate centroid of each slice
    other_dims = setdiff(1:3, iax);
    slice_idx = round(pts(:, iax));
    
    % Filter valid
    valid = slice_idx > 0 & slice_idx <= N;
    slice_idx = slice_idx(valid);
    pts = pts(valid, :);
    
    sum_vals = accumarray(slice_idx, 1, [N,1]);
    cx = accumarray(slice_idx, pts(:, other_dims(1)), [N,1]) ./ sum_vals;
    cy = accumarray(slice_idx, pts(:, other_dims(2)), [N,1]) ./ sum_vals;
    
    % Fill NaNs (empty slices)
    cx(sum_vals==0) = 0; cy(sum_vals==0) = 0;
    
    % The error is the distance from the global mean axis (0,0)
    % Assuming we want it centered at 0,0
    offsets = [cx, cy];
end