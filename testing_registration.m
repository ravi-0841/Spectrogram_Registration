warning('off', 'all');
[s_ang,~] = audioread('./angry.wav');
[s_neu,~] = audioread('./neutral.wav');

f = 16000;
r = 512;
w = 0.025;
o = 0.010;
linear_registration = 0;

[spect_ang_mag, spect_ang_phase] = get_spectrogram(s_ang,f,r,w,o);
[spect_neu_mag, spect_neu_phase] = get_spectrogram(s_neu,f,r,w,o);

n_cols = min([size(spect_neu_mag,2), size(spect_ang_mag,2)]);
res_spect_ang_mag = imresize(spect_ang_mag,[size(spect_ang_mag,1),n_cols]);
res_spect_neu_mag = imresize(spect_neu_mag,[size(spect_neu_mag,1),n_cols]);

res_spect_ang_phase = imresize(spect_ang_phase,[size(spect_ang_phase,1),n_cols]);
res_spect_neu_phase = imresize(spect_neu_phase,[size(spect_neu_phase,1),n_cols]);

%% Reconstruction Check
s_neu_recon = get_speech(spect_neu_mag,spect_neu_phase,f,r,w,o);
figure(), plot(s_neu, 'r'), hold on, plot(s_neu_recon, 'g'), title('Reconstruction');

%% Registration of magnitude spectrograms
if linear_registration          % Linear Registration
    figure(), subplot(121), imshowpair(res_spect_neu_mag,res_spect_ang_mag), title('Unregistered');
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.MaximumIterations = 300;
    transform = imregtform(res_spect_neu_mag,res_spect_ang_mag,'translation',optimizer,metric);
    registered_mag = imwarp(res_spect_neu_mag,transform);
    registered_phase = imwarp(res_spect_neu_phase,transform);
    subplot(122), imshowpair(registered_mag,res_spect_ang_mag), title('Linear Registration');
    res_reg_mag = imresize(registered_mag,size(spect_neu_mag));
    res_reg_phase = imresize(registered_phase,size(spect_neu_phase));
else                            % Non-Linear Registration
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.MaximumIterations = 300;
    transform = imregtform(res_spect_neu_mag,res_spect_ang_mag,'translation',optimizer,metric);
    res_spect_neu_mag = imwarp(res_spect_neu_mag,transform);
    res_spect_neu_phase = imwarp(res_spect_neu_phase,transform);
    
    recon_speech = get_speech(res_spect_neu_mag,res_spect_neu_phase,f,r,w,o);
    soundsc(recon_speech, f);
    
    [disp_field,movingReg] = imregdemons(res_spect_neu_mag,res_spect_ang_mag,[500,400,300],...
                                            'AccumulatedFieldSmoothing',2.5);
%     disp_field = block_demons(res_spect_neu_mag,res_spect_ang_mag,21,11,2.5);
    abs_disp_field = squeeze(disp_field(:,:,1)).^2 + squeeze(disp_field(:,:,2)).^2;
    
    registered_mag = imwarp(res_spect_neu_mag,disp_field);
    registered_phase = imwarp(res_spect_neu_phase,disp_field);
    
    figure(), subplot(121), imshowpair(res_spect_neu_phase,res_spect_ang_phase), title('Unregistered');
    subplot(122), imshowpair(registered_phase,res_spect_ang_phase), title('Phase Registration');
    
    figure(), subplot(121), imshowpair(res_spect_neu_mag,res_spect_ang_mag), title('Unregistered');
    subplot(122), imshowpair(registered_mag,res_spect_ang_mag), title('Spectrogram Registration');
    
    figure(), subplot(121), imshow(spect_neu_mag,[]), colormap('jet'), title('Original Spectrogram');
    subplot(122), imshow(registered_mag,[]), colormap('jet'), title('Modified Spectrogram');
    
    res_reg_mag = imresize(registered_mag,size(spect_neu_mag));
    res_reg_phase = imresize(registered_phase,size(spect_neu_phase));
    recon_speech = get_speech(res_reg_mag,res_reg_phase,f,r,w,o);
end

%% Visualize the vector field and apply independent transformations
if linear_registration==0
    flow = opticalFlow(-1*squeeze(disp_field(:,:,1)),-1*squeeze(disp_field(:,:,2)));
    figure(), subplot(121), plot(flow, 'DecimationFactor', [8,16]), title('Overall Displacement Field');
    subplot(122), imshow(abs_disp_field, []), colormap('jet'), title('Absolute Disp Field')

            % Warping only the frequency axis
    mod_disp_field = disp_field;
    mod_disp_field(:,:,1) = zeros(size(squeeze(mod_disp_field(:,:,1))));
    abs_warp = squeeze(mod_disp_field(:,:,1)).^2 + squeeze(mod_disp_field(:,:,2)).^2;

    registered_mag = imwarp(res_spect_neu_mag,mod_disp_field);
    registered_phase = imwarp(res_spect_neu_phase,mod_disp_field);

    flow = opticalFlow(-1*squeeze(mod_disp_field(:,:,1)),-1*squeeze(mod_disp_field(:,:,2)));
    figure(), subplot(121), plot(flow, 'DecimationFactor',[8,16]), title('Freq Axis Warping field');
    subplot(122), imshow(abs_warp, []), colormap('jet'), title('Absolute Warping field');

    figure(), subplot(121), imshowpair(res_spect_neu_phase,res_spect_ang_phase), title('Unregistered');
    subplot(122), imshowpair(registered_phase,res_spect_ang_phase), title('Freq Phase Registration');

    figure(), subplot(121), imshowpair(res_spect_neu_mag,res_spect_ang_mag), title('Unregistered');
    subplot(122), imshowpair(registered_mag,res_spect_ang_mag), title('Freq Spectrogram Registration');

    figure(), subplot(121), imshow(spect_neu_mag,[]), colormap('jet'), title('Original Spectrogram');
    subplot(122), imshow(registered_mag,[]), colormap('jet'), title('Freq Modified Spectrogram');

    res_reg_mag = imresize(registered_mag,size(spect_neu_mag));
    res_reg_phase = imresize(registered_phase,size(spect_neu_phase));
    recon_speech_freq_warped = get_speech(res_reg_mag,res_reg_phase,f,r,w,o);

            % Warping only the time axis
    mod_disp_field = disp_field;
    mod_disp_field(:,:,2) = zeros(size(squeeze(mod_disp_field(:,:,2))));
    abs_warp = squeeze(mod_disp_field(:,:,1)).^2 + squeeze(mod_disp_field(:,:,2)).^2;

    registered_mag = imwarp(res_spect_neu_mag,mod_disp_field);
    registered_phase = imwarp(res_spect_neu_phase,mod_disp_field);

    flow = opticalFlow(-1*squeeze(mod_disp_field(:,:,1)),-1*squeeze(mod_disp_field(:,:,2)));
    figure(), subplot(121), plot(flow, 'DecimationFactor',[8,16]), title('Time Axis Warping field');
    subplot(122), imshow(abs_warp, []), colormap('jet'), title('Absolute Warping field');

    figure(), subplot(121), imshowpair(res_spect_neu_phase,res_spect_ang_phase), title('Unregistered');
    subplot(122), imshowpair(registered_phase,res_spect_ang_phase), title('Time Phase Registration');

    figure(), subplot(121), imshowpair(res_spect_neu_mag,res_spect_ang_mag), title('Unregistered');
    subplot(122), imshowpair(registered_mag,res_spect_ang_mag), title('Time Spectrogram Registration');

    figure(), subplot(121), imshow(spect_neu_mag,[]), colormap('jet'), title('Original Spectrogram');
    subplot(122), imshow(registered_mag,[]), colormap('jet'), title('Time Modified Spectrogram');

    res_reg_mag = imresize(registered_mag,size(spect_neu_mag));
    res_reg_phase = imresize(registered_phase,size(spect_neu_phase));
    recon_speech_time_warped = get_speech(res_reg_mag,res_reg_phase,f,r,w,o);
end

%% Plot the original and reconstructed Speech
figure(), subplot(411), plot(s_neu), title('Original Speech')
subplot(412), plot(recon_speech), title('Reconstructed Speech')
subplot(413), plot(recon_speech_freq_warped), title('Freq warped Reconstructed Speech')
subplot(414), plot(recon_speech_time_warped), title('Time warped Reconstructed Speech')

%% Filtering with low pass filter
% fc = 4000;
% fs = 16000;
% [b,a] = butter(6,fc/(fs/2));
% filtered_speech = filter(b,a,recon_speech_sig);
% subplot(313), plot(filtered_speech), title('Filtered Reconstruction Speech')

%% Sound Playback
% soundsc(s_ang, f);
% pause(2)
% soundsc(recon_speech_sig, f);
% pause(2)
% soundsc(filtered_speech, f);


% close all





























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