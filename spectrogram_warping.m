% clear all
clc

N = 512;
wshift = 128;
Q=N/wshift;

W = hann(N);
S = hann(N);

%% For DTW
r = 512;
w = 0.025;
s = 0.010;
top = 3;

%% Get the wav files in
target = 'angry2.wav'; % please provide a test file  or Target
[x_tar,fs] = audioread(target);

source = 'neutral2.wav'; % please provide a test file  or Source
[x_src,fs] = audioread(source);

[x_src, x_tar] = get_alignment(x_src,x_tar,fs,w,w-s,r,top);

if size(x_src,2)>1, fprintf('This code only handles single channel files\n'); return; end
if size(x_tar,2)>1, fprintf('This code only handles single channel files\n'); return; end

fprintf(1,'Processing file: %s (length: %.2f s)\n',source,length(x_src)/fs);
fprintf(1,'Processing file: %s (length: %.2f s)\n',target,length(x_tar)/fs);

X_tar = stft(x_tar,N,wshift,W);
Xpow_tar = sum(abs(X_tar(:)).^2);

tmp_tar = stft(istft(abs(X_tar),wshift,S),N,wshift,W)-X_tar; % difference between temporary signal (zero phase) spectrogram and actual signal spectrogram 
C = 10*log10(Xpow_tar / sum(abs(tmp_tar(:)).^2));
fprintf(1,'Abs(X)          : %5.2f dB\n',C);

X_src = stft(x_src,N,wshift,W);
Xpow_src = sum(abs(X_src(:)).^2);

tmp_src = stft(istft(abs(X_src),wshift,S),N,wshift,W)-X_src; % difference between temporary signal (zero phase) spectrogram and actual signal spectrogram 
C = 10*log10(Xpow_src / sum(abs(tmp_src(:)).^2));
fprintf(1,'Abs(X)          : %5.2f dB\n',C);

%% Perform Demons Registration

X0_tar = abs(X_tar);
X0_src = abs(X_src);

opts = struct();
opts.alpha = 0.4;
opts.only_freq = 0;
opts.sigma_fluid = 1.5; %0.7
opts.sigma_diff = 2.5;  %1.5
opts.step = 1.0;
opts.max_iter = 600;
opts.pyramid_levels  = 1;
opts.compositive = 0;
opts.diffeomorphism = 1;
opts.plot = 1;

iterations = 1000;

disp_field = my_multires_demons(log(1+X0_tar),log(1+X0_src),opts);
warped_mag = imwarp(abs(X_src),disp_field);
warped_phase = imwarp(angle(X_src),disp_field);
recon_signal_fast = get_signal(warped_mag.*exp(1j*warped_phase),W,S,iterations,wshift);
recon_signal_iter = get_signal_iteratively(warped_mag.*exp(1j*warped_phase),N,wshift,W,iterations);

X_fast = stft(recon_signal_fast,N,wshift,W);
X_iter = stft(recon_signal_iter,N,wshift,W);

orig_phase = do_phase_unwrapping(angle(X_tar));
fast_phase = do_phase_unwrapping(angle(X_fast));
iter_phase = do_phase_unwrapping(angle(X_iter));

disp(['SSD orgVSfast: ' num2str(sum(sum((orig_phase - fast_phase).^2))) ...
    '    SSD orgVSiter: ' num2str(sum(sum((orig_phase - iter_phase).^2)))]);

figure()
lim = [1 1; size(warped_mag,1) size(warped_mag,2)];
subplot(131), imshowpair(log(1+X0_tar), log(1+warped_mag)), subplot(132), ...
    showgrid(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)),4,lim),...
    subplot(133), showvector(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)),5);

%% window wise registration along frequency axis
X0_tar = abs(X_tar);
X0_src = abs(X_src);

opts = struct();
opts.alpha = 0.4;
opts.only_freq = 0;
opts.sigma_fluid = 1.5; %0.7
opts.sigma_diff = 2.5;  %1.5
opts.step = 1.0;
opts.max_iter = 600;
opts.pyramid_levels  = 1;
opts.compositive = 0;
opts.diffeomorphism = 1;
opts.plot = 0;

iterations = 1000;

window_size = 7;
overlap = 3;
stride = window_size - overlap;
warped_spect = zeros(size(X0_tar));
for i = 1:stride:size(X0_tar,2)-window_size+1
    src_sub = X0_src(:, i:i+window_size-1);
    tar_sub = X0_tar(:, i:i+window_size-1);
    disp_field = my_multires_demons(log(1+tar_sub),log(1+src_sub),opts);
    warped_sub = imwarp(src_sub, disp_field);
    warped_spect(:, i:i+stride-1) = warped_sub(:,1:stride);
    if i == 1
        warped_spect(:, i:i+window_size-1) = warped_sub;
    else
        warped_spect(:, i+stride:i+window_size-1) = (warped_spect(:, ...
            i+stride:i+window_size-1) + warped_sub(:,stride+1:window_size))/2;
    end
end
recon_signal_fast = get_signal(warped_spect,W,S,iterations,wshift);