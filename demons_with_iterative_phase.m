clear all
clc

N=512;
wshift=128;
Q=N/wshift;

W = hann(N);
S = hann(N);

%% For DTW
r = 512;
w = 0.025;
s = 0.010;
top = 3;

%% Get the wav file in
infile1='angry.wav'; % please provide a test file
[x_ang,fs]=audioread(infile1);

infile2='neutral.wav'; % please provide a test file
[x_neu,fs]=audioread(infile2);

[x_neu, x_ang] = get_alignment(x_neu,x_ang,fs,w,w-s,r,top);

if size(x_ang,2)>1, fprintf('This code only handles single channel files\n'); return; end
if size(x_neu,2)>1, fprintf('This code only handles single channel files\n'); return; end

fprintf(1,'Processing file: %s (length: %.2f s)\n',infile1,length(x_ang)/fs);
fprintf(1,'Processing file: %s (length: %.2f s)\n',infile2,length(x_neu)/fs);

X_ang = stft(x_ang,N,wshift,W);
Xpow_ang = sum(abs(X_ang(:)).^2);

tmp_ang = stft(istft(abs(X_ang),wshift,S),N,wshift,W)-X_ang; % difference between temporary signal (zero phase) spectrogram and actual signal spectrogram 
C = 10*log10(Xpow_ang / sum(abs(tmp_ang(:)).^2));
fprintf(1,'Abs(X)          : %5.2f dB\n',C);

X_neu = stft(x_neu,N,wshift,W);
Xpow_neu = sum(abs(X_neu(:)).^2);

tmp_neu = stft(istft(abs(X_neu),wshift,S),N,wshift,W)-X_neu; % difference between temporary signal (zero phase) spectrogram and actual signal spectrogram 
C = 10*log10(Xpow_neu / sum(abs(tmp_neu(:)).^2));
fprintf(1,'Abs(X)          : %5.2f dB\n',C);

%% Perform Demons Registration

X0_ang = abs(X_ang);
X0_neu = abs(X_neu);

opts = struct();
opts.alpha = 0.4;
opts.sigma_fluid = 0.7;
opts.sigma_diff = 1.5;
opts.step = 1.0;
opts.compositive = 1;
opts.max_iter = 500;

do_diffeomorphic = 1;

if do_diffeomorphic
    disp_field = my_diffeomorphism(log(1 + X0_ang),log(1 + X0_neu),opts);
    warped_mag = imwarp(log(1 + abs(X_neu)),disp_field);
    warped_phase = imwarp(log(1 + angle(X_neu)),disp_field);
    est_signal_diffeo = get_signal_iteratively(warped_mag.*exp(1j*warped_phase), N, wshift, 1000);
    figure()
    lim = [1 1; size(warped_mag,1) size(warped_mag,2)];
    subplot(131), imshow(warped_mag,[]), colormap(jet), subplot(132), ...
        showgrid(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)),4,lim),...
        subplot(133), showvector(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)),5);
else
    disp_field = my_demons(log(1 + X0_ang),log(1 + X0_neu),opts);
    warped_mag = imwarp(log(1 + abs(X_neu)),disp_field);
    warped_phase = imwarp(angle(X_neu),disp_field);
    est_signal_demon = get_signal_iteratively(warped_mag.*exp(1j*warped_phase), N, wshift, 1000);
    figure()
    lim = [1 1; size(warped_mag,1) size(warped_mag,2)];
    subplot(131), imshow(warped_mag,[]), colormap(jet), subplot(132), ...
        showgrid(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)),4,lim), subplot(133), ...
        showvector(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)),5);
end