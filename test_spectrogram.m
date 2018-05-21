warning('off', 'all');
[s_ang,~] = audioread('./angry.wav');
[s_neu,~] = audioread('./neutral.wav');

f = 16000;
r = 512;
w = 0.025*f;
o = 0.010*f;
s = 0.015*f;

% [s_neu, s_ang] = get_alignment(s_neu, s_ang, f);

spect_ang = spectrogram(s_ang,w,o,r);
spect_neu = spectrogram(s_neu,w,o,r);

figure(), subplot(121), imshow(flipud(abs(spect_neu)), []), subplot(122), imshow(flipud(abs(spect_ang)), []), colormap('jet')

recon_ang_speech = get_speech(abs(spect_neu),angle(spect_neu),f,r,w/f,s/f,1);
soundsc(recon_ang_speech,f);
figure(), subplot(211), plot(s_neu), subplot(212), plot(recon_ang_speech)