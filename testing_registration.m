warning('off', 'all');
[s1,~] = audioread('./angry.wav');
[s2,~] = audioread('./neutral.wav');

f = 16000;
r = 512;
w = 0.025;
o = 0.010;
linear_registration = 1;

[spect_ang_mag, spect_ang_phase] = get_spectrogram(s1,f,r,w,o);
[spect_neu_mag, spect_neu_phase] = get_spectrogram(s2,f,r,w,o);

n_cols = min([size(spect_neu_mag,2), size(spect_ang_mag,2)]);
spect_ang = imresize(spect_ang_mag,[size(spect_ang_mag,1),n_cols]);
spect_neu = imresize(spect_neu_mag,[size(spect_neu_mag,1),n_cols]);

% spect_ang = imresize(spect_ang_phase,[size(spect_ang_phase,1),150]);
% spect_neu = imresize(spect_neu_phase,[size(spect_neu_phase,1),150]);

%% Registration of magnitude spectrograms
if linear_registration          % Linear transformations
    figure(), imshowpair(spect_neu,spect_ang), title('Unregistered')
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.MaximumIterations = 300;
    transform = imregtform(spect_neu,spect_ang,'translation',optimizer,metric);
    movingRegistered_mag = imwarp(spect_neu_mag,transform);
    movingRegistered_phase = imwarp(spect_neu_phase,transform);
    figure(), imshowpair(movingRegistered_mag,spect_ang), title('After Registration')
    transformed_mag = imresize(movingRegistered_mag,size(spect_neu_mag));
    transformed_phase = imresize(movingRegistered_phase,size(spect_neu_mag));
else             % Non-Linear Registration       
    figure(), imshowpair(spect_neu,spect_ang), title('Unregistered');
    [disp_field,movingReg] = imregdemons(spect_neu_mag,spect_ang_mag,100,...
                                            'AccumulatedFieldSmoothing',2);
    movingRegistered_mag = imwarp(spect_neu_mag,disp_field);
    movingRegistered_phase = imwarp(spect_neu_phase,disp_field);
    figure(), imshowpair(movingRegistered_mag,spect_ang), title('After Registration');
    transformed_mag = imresize(movingRegistered_mag,size(spect_neu_mag));
    transformed_phase = imresize(movingRegistered_phase,size(spect_neu_mag));
end

recon_speech_sig = get_speech(spect_neu_mag,spect_neu_phase,f,r,w,o);

%% Plot the original and reconstructed Speech
figure(), subplot(311), plot(s2), title('Original Speech')
subplot(312), plot(recon_speech_sig), title('Reconstructed Speech')

%% Filtering with low pass filter
fc = 4000;
fs = 16000;
[b,a] = butter(6,fc/(fs/2));
filtered_speech = filter(b,a,recon_speech_sig);
subplot(313), plot(filtered_speech), title('Filtered Reconstruction Speech')

%% Sound Playback
soundsc(s1, f);
pause(2)
soundsc(recon_speech_sig, f);
pause(2)
soundsc(filtered_speech, f);
































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