clear all

w = 0.025;
o = 0.010;
r = 512;

[d1,f] = audioread('./wav/angry.wav');
[d2,f] = audioread('./wav/neutral.wav');

d1 = -1 + 2*((d1 - min(d1)) / (max(d1) - min(d1)));
d2 = -1 + 2*((d2 - min(d2)) / (max(d2) - min(d2)));

ml = min(length(d1),length(d2));

% D1 = spectrogram(d1,w*f,int64(o*f),r);
% D2 = spectrogram(d2,w*f,int64(o*f),r);

D1 = specgram(d1,512,f,400,160);
D2 = specgram(d2,512,f,400,160);


[d1_feat,~,~] = mfcc(d1,f,25,15,0.97,'hamming',[130, 4000],20,26,22);
[d2_feat,~,~] = mfcc(d2,f,25,15,0.97,'hamming',[130, 4000],20,26,22);

SM = pair_sim(d1_feat,d2_feat,'cosine');
imagesc(1-SM), colormap(jet), colorbar;

[p,q,~] = dtw(1-SM, 'slopeThird');
hold on; plot(q,p,'w'); hold off

% Reconstruction using phase vocoder

% Calculate the frames in D2 that are indicated to match each frame in D1
idx_21 = zeros(1, size(D1,2));
for i = 1:length(idx_21) 
    idx_21(i) = q(find(p >= i,1)); 
end
% Phase-vocoder interpolate D2's STFT
D2x = pvsample(D2, idx_21-1, 128);
% Invert it back to time domain
d2x = istft(D2x, 512, 400, 240);
sound(d2x,f)

pause(2)

idx_12 = zeros(1, size(D2,2));
for i = 1:length(idx_12) 
    idx_12(i) = p(find(q >= i,1)); 
end
% Phase-vocoder interpolate D1's STFT
D1x = pvsample(D1, idx_12-1, 128);
% Invert it back to time domain
d1x = istft(D1x, 512, 400, 240);
sound(d1x,f)

%%
close all
clear all
addpath(genpath('./dtw'));
addpath(genpath('./mfcc'));

w = 0.025;
o = 0.010;
f = 16000;
r = 512;

[d1,~] = audioread('angry.wav');
[d2,~] = audioread('neutral.wav');

d1 = -1 + 2*((d1 - min(d1)) / (max(d1) - min(d1)));
d2 = -1 + 2*((d2 - min(d2)) / (max(d2) - min(d2)));

d1_spect = spectrogram(d1,int64(w*f),int64(o*f),r);
d2_spect = spectrogram(d2,int64(w*f),int64(o*f),r);

d1_mag = abs(d1_spect);
d1_phase = angle(d1_spect);

d2_mag = abs(d2_spect);
d2_phase = angle(d2_spect);

[d1_feat,~,~] = mfcc(d1,f,25,15,0.97,'hamming',[130, 4000],20,26,22);
[d2_feat,~,~] = mfcc(d2,f,25,15,0.97,'hamming',[130, 4000],20,26,22);

M = pair_sim(d1_feat, d2_feat, 'cosine');
figure(), imagesc(1-M), colormap(jet), colorbar;
[p,q,C] = dtw(1-M,'slopeThird');
hold on; plot(q,p,'w'); hold off

idx_21 = zeros(1, size(d1_mag,2));
for i = 1:length(idx_21) 
    idx_21(i) = q(find(p >= i,1)); 
end

idx_12 = zeros(1, size(d2_mag,2));
for i = 1:length(idx_12) 
    idx_12(i) = p(find(q >= i,1));
end

% Direct Reconstruction using stft interpolation
[d1_mag_tilda, d1_phase_tilda] = interpolate_stft_from_dtw(d1_mag,d1_phase,p); %idx_12
[d2_mag_tilda, d2_phase_tilda] = interpolate_stft_from_dtw(d2_mag,d2_phase,q); %idx_21

% figure(), subplot(121), imshow(d1_mag, []), title('D1 Original'), subplot(122), imshow(d1_mag_tilda, []), title('D1 Warped'), colormap('jet');
% figure(), subplot(121), imshow(d2_mag, []), title('D2 Original'), subplot(122), imshow(d2_mag_tilda, []), title('D2 Warped'), colormap('jet');
d1_recon = get_speech(d1_mag_tilda,d1_phase_tilda,f,512,0.025,0.015,1);
d2_recon = get_speech(d2_mag_tilda,d2_phase_tilda,f,512,0.025,0.015,1);
sound(d1_recon, f);
pause(2)
sound(d2_recon, f);