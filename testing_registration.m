warning('off', 'all');
[s1,~] = audioread('/home/ravi/Downloads/speech_data_male1/angry/51.wav');
[s2,~] = audioread('/home/ravi/Downloads/speech_data_male1/neutral/51.wav');

f = 16000;
r = 1024;
w = 0.025;
o = 0.010;

[spect_ang_org, spect_ang_phase] = get_spectrogram(s1,f,r,w,o);
[spect_neu_org, spect_neu_phase] = get_spectrogram(s2,f,r,w,o);

spect_ang = imresize(spect_ang_org,[size(spect_ang_org,1),100]);
spect_neu = imresize(spect_neu_org,[size(spect_neu_org,1),100]);

figure(), imshowpair(spect_neu,spect_ang), title('Unregistered'), colormap('jet')
[optimizer,metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 300;
movingRegisteredDefault = imregister(spect_neu,spect_ang,'affine',optimizer,metric);
figure(), imshowpair(movingRegisteredDefault,spect_ang), title('A: Default Registration')

transformed_image = imresize(movingRegisteredDefault,size(spect_neu_org));
recon_speech_sig = get_speech(transformed_image,spect_neu_phase,f,w,o);




%% Reconstruction using Sparsity
% A = dict;
% y = zeros(size(frames));
% recon_signal = zeros(length(s),1);
% for i=1:size(frames,2)
%     [reco_signal, ~] = omp(A,frames(:,i),'sparsity',9);
%     segment = A*reco_signal;
%     y(:,i) = segment;
%     if i==1
%         recon_sig(i:800*i) = segment;
%     else
%         recon_signal(400*(i-1)+1:400*i) = (recon_signal(400*(i-1)+1:400*i) + segment(1:400)) / 2;
%         recon_signal(400*i+1:400*(i+1)) = segment(401:end);
%     end
% end
% soundsc(recon_signal, 16000);