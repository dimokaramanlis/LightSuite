function [data, power_retained_spectrum] = spatial_bandpass_3d(data, radius, ...    
    f_lower_scale, f_upper_scale, use_gpu, smoothing_ratio)
% Butterworth bandpass filter for 3D images.

% Null op if filter cutoffs are trivial
if isinf(f_lower_scale) && isinf(f_upper_scale)
    power_retained_spectrum = 1;
    return;
end

% --- Smoothing Ratio Handling ---
if ~exist('smoothing_ratio', 'var')
    scale_y = 1;
    scale_x = 1;
    scale_z = 1;
elseif numel(smoothing_ratio)==1
     if smoothing_ratio > 1
        scale_x = smoothing_ratio;
        scale_y = 1;
        scale_z = 1;
     else
        scale_x = 1;
        scale_y = 1/smoothing_ratio;
        scale_z = 1;
     end
elseif numel(smoothing_ratio) == 3
    scale_x = smoothing_ratio(1);
    scale_y = smoothing_ratio(2);
    scale_z = smoothing_ratio(3);
else
    error('smoothing_ratio must be a scalar or a 3-element vector');
end

% --- GPU Handling ---
if ~exist('use_gpu', 'var')
    use_gpu = isa(data, 'gpuArray');
end
is_input_gpuArray = isa(data, 'gpuArray');

% --- Data Dimensions ---
[h, w, t] = size(data);

% --- FFT Dimensions (Corrected) ---
hf = single(2^nextpow2(h));  % FFT size for height
wf = single(2^nextpow2(w));  % FFT size for width
tf = single(2^nextpow2(t));  % FFT size for time/depth
if use_gpu
    hf = gpuArray(hf);  % FFT size for height
    wf = gpuArray(wf);  % FFT size for width
    tf = gpuArray(tf);  % FFT size for time/depth
end

% --- Butterworth Filter Parameters ---
n = 4; % degree of the Butterworth polynomial
sigma = radius/2;
f_corner = 1 / pi ./ sigma;
f_lower = f_corner / f_lower_scale;
f_upper = f_corner * f_upper_scale;

% --- Create Butterworth Filter ---
[cx, cy, cz] = meshgrid(1:wf, 1:hf, 1:tf);
dist_matrix = sqrt( (scale_y*((cy-1) / (hf) - 1 / 2)).^2 + ...
                    (scale_x*((cx-1) / (wf) - 1 / 2)).^2 + ...
                    (scale_z*((cz-1) / (tf) - 1 / 2)).^2);

dist_matrix = max(dist_matrix, 1e-6); % Avoid NaN

lpf = 1 ./ (1 + (dist_matrix / f_upper).^(2 * n));
hpf = 1 - 1 ./ (1 + (dist_matrix / f_lower).^(2 * n));
bpf = single(lpf .* hpf);  % Convert to single precision
clear lpf hpf dist_matrix cx cy cz
power_retained_spectrum = (sum(bpf(:).^2) / (hf*wf*tf));
bpf = maybe_gpu(use_gpu, bpf);

% --- Memory Management (Robust) ---

% Calculate the size of the data and filter (in GB)
data_size_gb = whos('data').bytes / (1024^3);
bpf_size_gb = whos('bpf').bytes / (1024^3);  % Even if on GPU, whos() works
total_memory_needed_gb = data_size_gb + bpf_size_gb + data_size_gb; % + output

% Get available memory (in GB)
available_memory_gb = get_free_mem() / (1024^3);

% Check if enough memory is available
if use_gpu
    gpu_info = gpuDevice();
    available_memory_gb = gpu_info.AvailableMemory / (1024^3);
    if total_memory_needed_gb * 2 > available_memory_gb  % x2 safety factor
      error(['Not enough GPU memory.  Need approx. ', num2str(total_memory_needed_gb*2), ' GB, but only ', num2str(available_memory_gb), ' GB is available.  Consider reducing data size or using CPU (use_gpu = false).']);
    end
else
    if total_memory_needed_gb * 4 > available_memory_gb   % x4 safety factor for CPU
        error(['Not enough RAM. Need approx. ', num2str(total_memory_needed_gb*4), ' GB, but only ', num2str(available_memory_gb), ' GB is available. Consider reducing data size or increasing RAM.']);
    end
end


% --- Main Processing ---
data = maybe_gpu(use_gpu & ~is_input_gpuArray, data);

% pad
pad_each_h = floor((hf-h)/2);
pad_each_w = floor((wf-w)/2);
pad_each_t = floor((tf-t)/2);
data = padarray(data, double([pad_each_h, pad_each_w, pad_each_t]),'symmetric', 'both');

% fft
data = fftn(data, [hf, wf, tf]);
data = fftshift(data);

% Filter
data = data .* bpf;

data = ifftshift(data);

% ifft
data = ifftn(data);
data = real(data);

% unpad
data = data((pad_each_h+1):(pad_each_h+h), ...
            (pad_each_w+1):(pad_each_w+w), ...
            (pad_each_t+1):(pad_each_t+t));

% Move back to CPU
if use_gpu && ~is_input_gpuArray
    data = gather(data);
end

end

% --- Helper Functions ---
function x = maybe_gpu(use_gpu, x)
    if use_gpu
        x = gpuArray(x);
    end
end

function f = get_free_mem()
% Get free memory in bytes (cross-platform)
if ispc
    [~, userview] = memory;
    f = userview.PhysicalMemory.Available;
elseif ismac
    [~, mem_stats] = unix('vm_stat');
    mem_stats = strsplit(mem_stats, '\n');
    mem_stats = cellfun(@(s) str2double(regexp(s, '[0-9]+', 'match')), mem_stats(2:5), 'UniformOutput', false);
     if ~isempty(mem_stats) && ~any(cellfun(@isempty, mem_stats))
        mem_stats = cell2mat(mem_stats);
        f = mem_stats(1) * mem_stats(4) * 4096;
     else
         f = 4e9; % Return default 4GB
     end
else
    [~, w] = unix('free -b');
    stats = strsplit(w, '\n');
    stats = strsplit(stats{2}, ' ');
    f = str2double(stats{end-1});
end
end
