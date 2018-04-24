[s1,~] = audioread('/home/ravi/Downloads/speech_data_male1/angry/50.wav');
[s2,~] = audioread('/home/ravi/Downloads/speech_data_male1/neutral/50.wav');

f = 16000;
w = 0.025;
o = 0.010;

[spect_ang_org, spect_ang_phase] = get_spectrogram(s1,f,w,o);
[spect_neu_org, spect_neu_phase] = get_spectrogram(s2,f,w,o);

spect_ang = imresize(spect_ang_org,[size(spect_ang,1),100]);
spect_neu = imresize(spect_neu_org,[size(spect_neu,1),100]);

figure(), imshowpair(spect_ang,spect_neu), title('Unregistered'), colormap('jet')
[optimizer,metric] = imregconfig('multimodal');
movingRegisteredDefault = imregister(spect_neu,spect_ang,'affine',optimizer,metric);
figure(), imshowpair(movingRegisteredDefault,spect_ang), title('A: Default Registration')

transformed_image = imresize(movingRegisteredDefault,size(spect_neu_org));
recon_speech_sig = get_speech(transformed_image,spect_neu_phase,f,w,o);