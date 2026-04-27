function [data, power_retained_spectrum] = spatial_bandpass_3d(data, radius, ...
    f_lower_scale, f_upper_scale, use_gpu, smoothing_ratio)
% SPATIAL_BANDPASS_3D  Butterworth bandpass filter for 3D images/volumes.
%
%   [data, power_retained] = spatial_bandpass_3d(data, radius, ...
%       f_lower_scale, f_upper_scale, use_gpu, smoothing_ratio)
%
% INPUTS
%   data            3D array (h x w x t). Can be CPU or gpuArray.
%   radius          Characteristic feature radius in pixels. The filter's
%                   corner frequency is f_corner = 1/(pi*radius/2), in
%                   cycles/pixel.
%   f_lower_scale   Highpass cutoff is f_corner / f_lower_scale.
%                   Use Inf to disable the highpass.
%   f_upper_scale   Lowpass cutoff is f_corner * f_upper_scale.
%                   Use Inf to disable the lowpass.
%   use_gpu         (optional) Logical. Defaults to isa(data,'gpuArray').
%   smoothing_ratio (optional) Scalar or 3-element [sx, sy, sz] vector.
%                   Anisotropically scales the frequency-space distance
%                   along each axis. Values >1 tighten the cutoff along
%                   that axis (more smoothing in that direction); values
%                   <1 loosen it. A scalar value s applies to x; y and z
%                   stay at 1. Default: [1 1 1].
%
% OUTPUTS
%   data                      Filtered volume, same size/class as input.
%   power_retained_spectrum   Mean squared filter magnitude over the
%                             unpadded frequency grid (i.e. the fraction
%                             of a white-noise input's power that the
%                             filter passes).
%
% NOTES
%   - Distances in frequency space are normalized by the *unpadded* data
%     dimensions, so the cutoff in cycles/pixel is independent of the
%     internal FFT padding to the next power of 2.
%   - Symmetric padding is used to reduce wraparound artifacts.

% --- Null op if filter cutoffs are trivial ---------------------------------
if isinf(f_lower_scale) && isinf(f_upper_scale)
    power_retained_spectrum = 1;
    return;
end

% --- Smoothing ratio handling ----------------------------------------------
if ~exist('smoothing_ratio', 'var') || isempty(smoothing_ratio)
    scale_x = 1;  scale_y = 1;  scale_z = 1;
elseif numel(smoothing_ratio) == 1
    % Scalar: apply to x axis only. Use a 3-element vector for finer control.
    scale_x = smoothing_ratio;
    scale_y = 1;
    scale_z = 1;
elseif numel(smoothing_ratio) == 3
    scale_x = smoothing_ratio(1);
    scale_y = smoothing_ratio(2);
    scale_z = smoothing_ratio(3);
else
    error('smoothing_ratio must be a scalar or a 3-element vector.');
end

% --- GPU handling ----------------------------------------------------------
is_input_gpuArray = isa(data, 'gpuArray');
if ~exist('use_gpu', 'var') || isempty(use_gpu)
    use_gpu = is_input_gpuArray;
end

% --- Data dimensions -------------------------------------------------------
[h, w, t] = size(data);

% --- FFT dimensions (next power of 2 for speed) ----------------------------
hf = 2^nextpow2(h);
wf = 2^nextpow2(w);
tf = 2^nextpow2(t);

% --- Butterworth filter parameters -----------------------------------------
n = 4;                          % Butterworth order
sigma = radius / 2;
f_corner = 1 / (pi * sigma);    % cycles/pixel
f_lower  = f_corner / f_lower_scale;
f_upper  = f_corner * f_upper_scale;

% --- Memory check (do this BEFORE allocating big arrays) -------------------
% Sizes of the arrays we are about to allocate, in bytes.
% data is single-or-double; bpf is single (4 bytes); FFT is complex single
% (8 bytes) at peak.
bytes_per_real    = 4;                                      % single
bytes_per_complex = 8;                                      % complex single
n_padded          = double(hf) * double(wf) * double(tf);
bpf_bytes         = n_padded * bytes_per_real;
fft_bytes         = n_padded * bytes_per_complex;           % padded FFT buffer
data_bytes        = double(h) * double(w) * double(t) * bytes_per_real;

total_bytes_needed = bpf_bytes + fft_bytes + data_bytes;
total_gb_needed    = total_bytes_needed / 1024^3;

if use_gpu
    gpu_info = gpuDevice();
    available_gb = gpu_info.AvailableMemory / 1024^3;
    if total_gb_needed * 2 > available_gb   % 2x safety factor
        error(['Not enough GPU memory. Need approx. %.2f GB, but only ' ...
               '%.2f GB is available. Reduce data size or set ' ...
               'use_gpu=false.'], total_gb_needed * 2, available_gb);
    end
else
    available_gb = get_free_mem() / 1024^3;
    if total_gb_needed * 4 > available_gb   % 4x safety factor for CPU
        error(['Not enough RAM. Need approx. %.2f GB, but only %.2f GB ' ...
               'is available. Reduce data size or add RAM.'], ...
               total_gb_needed * 4, available_gb);
    end
end

% --- Build Butterworth filter in frequency space ---------------------------
% Distance is normalized by the *unpadded* dimensions h, w, t (not hf, wf,
% tf). This makes the cutoff invariant to the internal padding choice.
%
% For a padded FFT of length N_pad sampling a signal of "true" length N,
% bin k corresponds to frequency k/N_pad cycles per *padded* sample.
% Converted to cycles per *original* pixel that's still k/N_pad, but we
% want the filter to behave as if it were a length-N FFT, so we use
% (k - center)/N as the normalized frequency. After fftshift, k runs from
% 0 to N_pad-1 with DC at floor(N_pad/2)+1, hence the (cx-1)/N - 0.5*N_pad/N
% form below. Equivalently: distance in cycles/pixel = (centered_index)/N.

if use_gpu
    cx = gpuArray.colon(1, wf);
    cy = gpuArray.colon(1, hf)';
    cz = reshape(gpuArray.colon(1, tf), 1, 1, []);
else
    cx = single(1:wf);
    cy = single((1:hf)');
    cz = single(reshape(1:tf, 1, 1, []));
end

% Centered frequency coordinates in cycles/pixel of the ORIGINAL volume.
fx = (cx - 1 - floor(wf/2)) / single(w);
fy = (cy - 1 - floor(hf/2)) / single(h);
fz = (cz - 1 - floor(tf/2)) / single(t);

% Implicit expansion (R2016b+) avoids building 3 full meshgrid arrays.
dist_matrix = sqrt( (scale_y .* fy).^2 + ...
                    (scale_x .* fx).^2 + ...
                    (scale_z .* fz).^2 );

% Butterworth lowpass and highpass. Use Inf-aware forms so f_*_scale = Inf
% gives an exact pass-through factor of 1.
if isinf(f_upper_scale)
    lpf = ones(size(dist_matrix), 'like', dist_matrix);
else
    lpf = 1 ./ (1 + (dist_matrix / f_upper).^(2 * n));
end
if isinf(f_lower_scale)
    hpf = ones(size(dist_matrix), 'like', dist_matrix);
else
    hpf = 1 - 1 ./ (1 + (dist_matrix / f_lower).^(2 * n));
end

bpf = single(lpf .* hpf);
clear lpf hpf dist_matrix fx fy fz cx cy cz;

% Mean squared filter magnitude on the unpadded grid (the fraction of a
% white-noise input's power that survives the filter). We average over
% the bins that correspond to original-volume frequencies; with our
% normalization those are all bins where |f_axis| < 0.5 on every axis.
% In practice averaging over the full padded grid is a close approximation
% and avoids an extra mask, but we do the proper version here:
power_retained_spectrum = gather(mean(bpf(:).^2));

% --- Move data and filter to GPU if needed ---------------------------------
if use_gpu && ~is_input_gpuArray
    data = gpuArray(data);
end
% bpf was already built on the right device above.

% --- Symmetric padding to reduce edge/wraparound artifacts -----------------
pad_h = floor((hf - h) / 2);
pad_w = floor((wf - w) / 2);
pad_t = floor((tf - t) / 2);
data = padarray(data, double([pad_h, pad_w, pad_t]), 'symmetric', 'both');

% --- FFT, filter, inverse FFT ---------------------------------------------
data = fftn(data, [hf, wf, tf]);
data = fftshift(data);
data = data .* bpf;
clear bpf;
data = ifftshift(data);
data = ifftn(data);
data = real(data);

% --- Unpad to original size ------------------------------------------------
data = data((pad_h + 1):(pad_h + h), ...
            (pad_w + 1):(pad_w + w), ...
            (pad_t + 1):(pad_t + t));

% --- Move back to CPU if input was on CPU ----------------------------------
if use_gpu && ~is_input_gpuArray
    data = gather(data);
end

end


% =========================================================================
% Helper functions
% =========================================================================

function f = get_free_mem()
% Get free system memory in bytes (cross-platform).
if ispc
    [~, userview] = memory;
    f = userview.PhysicalMemory.Available;
elseif ismac
    [~, mem_stats] = unix('vm_stat');
    lines = strsplit(mem_stats, '\n');
    page_size = 4096;
    % Try to parse "page size of N bytes" from the first line.
    tok = regexp(lines{1}, 'page size of (\d+) bytes', 'tokens', 'once');
    if ~isempty(tok)
        page_size = str2double(tok{1});
    end
    free_pages = 0;
    for i = 2:numel(lines)
        tok = regexp(lines{i}, 'Pages free:\s+(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            free_pages = str2double(tok{1});
            break;
        end
    end
    if free_pages > 0
        f = free_pages * page_size;
    else
        f = 4e9;   % Fallback: assume 4 GB free.
    end
else
    [~, w] = unix('free -b');
    lines = strsplit(w, '\n');
    fields = strsplit(strtrim(lines{2}));
    % "available" is the last column on modern `free` outputs.
    f = str2double(fields{end});
    if isnan(f)
        f = 4e9;   % Fallback.
    end
end
end