%% Pipeline
do_no_future_lws = 1;
do_online_lws    = 1;
do_batch_lws     = 1;

%% Wav output flag
output_wav = 1;
%% Set analysis parameters
N=512;
wshift=128;
Q=N/wshift;
% Note: the code assumes the windows are symmetric, i.e., W(n) = W(N+1-n)
%W=sqrt((0.5-0.5*cos(2*pi*(1:2:2*N-1)'/(2*N)))/Q*2);
%S=sqrt((0.5-0.5*cos(2*pi*(1:2:2*N-1)'/(2*N)))/Q*2);

% W=sqrt((0.5-0.5*cos(2*pi*(0:(N-1))'/(N)))/Q*2);
% S=sqrt((0.5-0.5*cos(2*pi*(0:(N-1))'/(N)))/Q*2);

W = hann(N);
S = hann(N);
%% Set DTW Parameters
r = 512;
w = 0.025;
s = 0.010;
top = 3;

%% Get the wav file in
infile1='neutral.wav'; % please provide a test file
[x_ang,fs]=audioread(infile1);

infile2='angry.wav'; % please provide a test file
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

%% Create weights
L=5;
weights=create_weights(W,S,wshift,L);

%% Create weights for asymmetric windows
[W_asym_init,W_asym_full] = build_asymmetric_windows(W.*S,wshift);
weights_asym_init=create_weights(W_asym_init,S,wshift,L);
weights_asym_full=create_weights(W_asym_full,S,wshift,L);

%% Perform Demons Registration
% At this point, this does not seem to help much w.r.t. doing nothing

X0_ang = abs(X_ang);
X0_neu = abs(X_neu);

[disp_field, ~] = my_demons(log(1 + X0_ang),log(1 + X0_neu),0.4,0.7,1.5,1,0.000001);
% disp_field(:,:,1) = zeros(size(squeeze(disp_field(:,:,1))));
warped_mag = imwarp(log(1 + abs(X_neu)),disp_field);
warped_phase = imwarp(angle(X_neu),disp_field);

%% Perform a few iterations of Online LWS (aka. TF-RTISI-LA)
if do_online_lws
    tic;
    look_ahead        = 3;
    online_iterations = 100;
    online_alpha      = 1;
    online_beta       = 0.1;
    online_gamma      = 1;
    online_thresholds = online_alpha*exp(-online_beta*(0:(online_iterations-1)).^online_gamma);
    
    X1=online_lws(warped_mag.*exp(1j*angle(warped_phase)),weights,weights_asym_init,weights_asym_full,online_thresholds,look_ahead);
    time_online_lws = toc;
    
    tmp = stft(istft(X1,wshift,S),N,wshift,W)-X1;
    C1 = 10*log10(Xpow_neu / sum(abs(tmp(:)).^2));
    fprintf(1,'+Online LWS     : %5.2f dB (time: %.2f s)\n',C1,time_online_lws);
    if output_wav
        x1=istft(X1,wshift,S);
        audiowrite('out_online_lws.wav',x1,fs);
    end
else
    X1 = X0_neu;
end

%% run batch LWS
if do_batch_lws
    tic;
    iterations = 1000;
    alpha      = 100;
    beta       = 0.1;
    gamma      = 1;
    thresholds = alpha*exp(-beta*(0:(iterations-1)).^gamma);

    Y=batch_lws(X1,weights,thresholds);
    time_bach_lws = toc;
    
    tmp = stft(istft(Y,wshift,S),N,wshift,W)-Y;
    Cfinal = 10*log10(Xpow_neu / sum(abs(tmp(:)).^2));
    fprintf(1,'+Batch LWS      : %5.2f dB (time: %.2f s)\n',Cfinal,time_bach_lws);
    if output_wav
        y=istft(Y,wshift,S);
        audiowrite('out_batch_lws.wav',y,fs);
    end

else
    Y = X1;
end