[d1,sr] = audioread('angry.wav');
[d2,sr] = audioread('neutral.wav');

d1 = -1 + 2*((d1 - min(d1)) / (max(d1) - min(d1)));
d2 = -1 + 2*((d2 - min(d2)) / (max(d2) - min(d2)));
 
% Listen to them together:
ml = min(length(d1),length(d2));
% soundsc(d1(1:ml)+d2(1:ml),sr)
% or, in stereo
% soundsc([d1(1:ml),d2(1:ml)],sr)

% Calculate STFT features for both sounds (25% window overlap)
D1 = specgram(d1,512,sr,400,240);
D2 = specgram(d2,512,sr,400,240);

% Construct the 'local match' scores matrix as the cosine distance 
% between the STFT magnitudes
SM = simmx(abs(D1),abs(D2));
% Look at it:
subplot(121)
imagesc(SM)
colormap(1-gray)
% You can see a dark stripe (high similarity values) approximately
% down the leading diagonal.

% Use dynamic programming to find the lowest-cost path between the 
% opposite corners of the cost matrix
% Note that we use 1-SM because dp will find the *lowest* total cost
[p,q,C] = dp(1-SM);
% Overlay the path on the local similarity matrix
hold on; plot(q,p,'r'); hold off
% Path visibly follows the dark stripe

% Plot the minimum-cost-to-this point matrix too
subplot(122)
imagesc(C)
hold on; plot(q,p,'r'); hold off

% Bottom right corner of C gives cost of minimum-cost alignment of the two
C(size(C,1),size(C,2))
% This is the value we would compare between different 
% templates if we were doing classification.

% Calculate the frames in D2 that are indicated to match each frame
% in D1, so we can resynthesize a warped, aligned version
D2i1 = zeros(1, size(D1,2));
for i = 1:length(D2i1); D2i1(i) = q(find(p >= i,1)); end
% Phase-vocoder interpolate D2's STFT under the time warp
D2x = pvsample(D2, D2i1-1, 128);
% Invert it back to time domain
d2x = istft(D2x, 512, 512, 128);

% Listen to the results
% Warped version alone
sound(d2x,sr)
% Warped version added to original target (have to fine-tune length)
% d2x = resize(d2x', length(d1),1);
% soundsc(d1+d2x,sr)
% .. and in stereo
% soundsc([d1,d2x],sr)
% Compare to unwarped pair:
% soundsc([d1(1:ml),d2(1:ml)],sr)

%%
close all
clear all
addpath(genpath('./dtw'));
addpath(genpath('./mfcc'));

[d1,sr] = audioread('angry.wav');
[d2,sr] = audioread('neutral.wav');

d1 = -1 + 2*((d1 - min(d1)) / (max(d1) - min(d1)));
d2 = -1 + 2*((d2 - min(d2)) / (max(d2) - min(d2)));

[d1_mag, d1_phase] = get_spectrogram(d1,sr,512,0.025,0.015,0);
[d2_mag, d2_phase] = get_spectrogram(d2,sr,512,0.025,0.015,0);

[d1_feat,~,~] = mfcc(d1,sr,25,15,0.97,'hamming',[130, 4000],20,26,22);
[d2_feat,~,~] = mfcc(d2,sr,25,15,0.97,'hamming',[130, 4000],20,26,22);

M = simmx(d1_feat, d2_feat, 'cosine');
figure(), imshow(M, [])
[p,q,C] = dp(1-M);
hold on; plot(q,p,'r'); hold off

d1_mag_tilda = d1_mag(:,p);
d1_phase_tilda = d1_phase(:,p);

d2_mag_tilda = d2_mag(:,q);
d2_phase_tilda = d2_phase(:,q);

figure(), subplot(121), imshow(d1_mag, []), title('D1 Original'), subplot(122), imshow(d1_mag_tilda, []), title('D1 Warped'), colormap('jet');
figure(), subplot(121), imshow(d2_mag, []), title('D2 Original'), subplot(122), imshow(d2_mag_tilda, []), title('D2 Warped'), colormap('jet');
d1_recon = get_speech(d1_mag_tilda,d1_phase_tilda,sr,512,0.025,0.015);
d2_recon = get_speech(d2_mag_tilda,d2_phase_tilda,sr,512,0.025,0.015);
sound(d1_recon, sr);
pause(2)
sound(d2_recon, sr);