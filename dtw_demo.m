clear all

w = 0.025;
o = 0.015;
r = 512;

[d1,f] = audioread('angry.wav');
[d2,f] = audioread('neutral.wav');

d1 = -1 + 2*((d1 - min(d1)) / (max(d1) - min(d1)));
d2 = -1 + 2*((d2 - min(d2)) / (max(d2) - min(d2)));

ml = min(length(d1),length(d2));

D1 = spectrogram(d1,w*f,int64(o*f),r);
D2 = spectrogram(d2,w*f,int64(o*f),r);

d1_mag = abs(D1);
d1_phase = angle(D1);

d2_mag = abs(D2);
d2_phase = angle(D2);

[d1_feat,~,~] = mfcc(d1,f,25,10,0.97,'hamming',[130, 4000],20,26,22);
[d2_feat,~,~] = mfcc(d2,f,25,10,0.97,'hamming',[130, 4000],20,26,22);

SM = pair_sim(d1_feat,d2_feat,'cosine');
subplot(121)
imagesc(SM)
colormap(1-gray)

[p,q,C] = dtw(1-SM, 'slopeThird');
hold on; plot(q,p,'r'); hold off
subplot(122)
imagesc(C)
hold on; plot(q,p,'r'); hold off

% Bottom right corner of C gives cost of minimum-cost alignment of the two
C(size(C,1),size(C,2))

% Calculate the frames in D2 that are indicated to match each frame
% in D1, so we can resynthesize a warped, aligned version
D2i1 = zeros(1, size(D1,2));
for i = 1:length(D2i1) 
    D2i1(i) = q(find(p >= i,1)); 
end
% Phase-vocoder interpolate D2's STFT under the time warp
D2x = pvsample(D2, D2i1-1, 128);
% Invert it back to time domain
d2x = istft(D2x, 512, 512, 128);
sound(d2x,f)

pause(2)

D1i1 = zeros(1, size(D2,2));
for i = 1:length(D1i1) 
    D1i1(i) = p(find(q >= i,1)); 
end
% Phase-vocoder interpolate D2's STFT under the time warp
D1x = pvsample(D1, D1i1-1, 128);
% Invert it back to time domain
d1x = istft(D1x, 512, 512, 128);
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

d1_spect = spectrogram(d1,w*f,int64(o*f),r);
d2_spect = spectrogram(d2,w*f,int64(o*f),r);

d1_mag = abs(d1_spect);
d1_phase = angle(d1_spect);

d2_mag = abs(d2_spect);
d2_phase = angle(d2_spect);

[d1_feat,~,~] = mfcc(d1,f,25,15,0.97,'hamming',[130, 4000],20,26,22);
[d2_feat,~,~] = mfcc(d2,f,25,15,0.97,'hamming',[130, 4000],20,26,22);

M = pair_sim(d1_feat, d2_feat, 'cosine');
figure(), imshow(M, [])
[p,q,C] = dtw(1-M,'slopeThird');
hold on; plot(q,p,'r'); hold off

% d1_mag_tilda = d1_mag(:,p);
% d1_phase_tilda = d1_phase(:,p);
% 
% d2_mag_tilda = d2_mag(:,q);
% d2_phase_tilda = d2_phase(:,q);

[d1_mag_tilda, d1_phase_tilda] = interpolate_stft_from_dtw(d1_mag,d1_phase,p);
[d2_mag_tilda, d2_phase_tilda] = interpolate_stft_from_dtw(d2_mag,d2_phase,q);

figure(), subplot(121), imshow(d1_mag, []), title('D1 Original'), subplot(122), imshow(d1_mag_tilda, []), title('D1 Warped'), colormap('jet');
figure(), subplot(121), imshow(d2_mag, []), title('D2 Original'), subplot(122), imshow(d2_mag_tilda, []), title('D2 Warped'), colormap('jet');
d1_recon = get_speech(d1_mag_tilda,d1_phase_tilda,f,512,0.025,0.015,1);
d2_recon = get_speech(d2_mag_tilda,d2_phase_tilda,f,512,0.025,0.015,1);
sound(d1_recon, f);
pause(2)
sound(d2_recon, f);