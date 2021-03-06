warning('off', 'all');
[s_ang,~] = audioread('./wav/angry.wav');
[s_neu,~] = audioread('./wav/neutral.wav');

f = 16000;
r = 512;
w = 0.025;
s = 0.010;
linear_registration = 0;

[s_neu, s_ang] = get_alignment(s_neu,s_ang,f,w,w-s,r,3);

spect_ang = spectrogram(s_ang,w*f,int64((w-s)*f),r);
spect_neu = spectrogram(s_neu,w*f,int64((w-s)*f),r);

% spect_ang = imresize(spect_ang, [r/2 + 1, min([size(spect_ang,2), size(spect_neu,2)])]);
% spect_neu = imresize(spect_neu, [r/2 + 1, min([size(spect_ang,2), size(spect_neu,2)])]);

spect_ang_mag = abs(spect_ang);
spect_ang_phase = angle(spect_ang);

spect_neu_mag = abs(spect_neu);
spect_neu_phase = angle(spect_neu);

spect_neu_mag = mat2gray(log(1 + spect_neu_mag));
spect_ang_mag = mat2gray(log(1 + spect_ang_mag));

figure(), subplot(121), imshow(spect_neu_mag, []), title('Neutral'), subplot(122), imshow(spect_ang_mag, []), title('Angry'), colormap('jet')

%% Reconstruction Check
s_ang_recon = get_speech(spect_ang_mag,spect_ang_phase,f,r,w,s,1);
s_ang_recon = -1 + 2*(s_ang_recon - min(s_ang_recon))/(max(s_ang_recon) - min(s_ang_recon));
figure(1), plot(s_ang, 'r')
figure(1), hold on, plot(s_ang_recon, 'g'), title('Reconstruction');

%% Registration of magnitude spectrograms
if linear_registration          % Linear Registration
    figure(), subplot(121), imshowpair(res_spect_neu_mag,res_spect_ang_mag), title('Unregistered');
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.MaximumIterations = 300;
    transform = imregtform(normc(spect_neu_mag),normc(spect_ang_mag),'translation',optimizer,metric);
    registered_mag = imwarp(spect_neu_mag,transform);
    registered_phase = imwarp(spect_neu_phase,transform);
    subplot(122), imshowpair(registered_mag,spect_ang_mag), title('Linear Registration');
%     res_reg_mag = imresize(registered_mag,size(spect_neu_mag));
%     res_reg_phase = imresize(registered_phase,size(spect_neu_phase));

else                            % Non-Linear Registration
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.MaximumIterations = 300;
    transform = imregtform(normc(spect_neu_mag),normc(spect_ang_mag),'translation',optimizer,metric);
    spect_neu_mag = imwarp(spect_neu_mag,transform);
    spect_neu_phase = imwarp(spect_neu_phase,transform);
    
    recon_speech = get_speech(spect_neu_mag,spect_neu_phase,f,r,w,s,1);
    soundsc(recon_speech, f);
    
    [disp_field,movingReg] = imregdemons(normc(spect_neu_mag),normc(spect_ang_mag),[500,400,300],...
                                            'AccumulatedFieldSmoothing',0.005,'PyramidLevels',3,...
                                            'DisplayWaitbar', false);
                                        
%     disp_field = block_demons(res_spect_neu_mag,res_spect_ang_mag,21,11,2.5);
    abs_disp_field = squeeze(disp_field(:,:,1)).^2 + squeeze(disp_field(:,:,2)).^2;
    
    registered_mag = imwarp(spect_neu_mag,disp_field);
    registered_phase = imwarp(spect_neu_phase,disp_field);
    
%     figure(), subplot(121), imshowpair(res_spect_neu_phase,res_spect_ang_phase), title('Unregistered');
%     subplot(122), imshowpair(registered_phase,res_spect_ang_phase), title('Phase Registration');
    
    figure(), subplot(121), imshowpair(spect_neu_mag,spect_ang_mag), title('Unregistered');
    subplot(122), imshowpair(registered_mag,spect_ang_mag), title('Spectrogram Registration');
    
    figure(), subplot(121), imshow(spect_neu_mag,[]), colormap('jet'), title('Original Spectrogram');
    subplot(122), imshow(registered_mag,[]), colormap('jet'), title('Modified Spectrogram');
    
%     res_reg_mag = imresize(registered_mag,size(spect_neu_mag));
%     res_reg_phase = imresize(registered_phase,size(spect_neu_phase));
    recon_speech = get_speech(registered_mag,registered_phase,f,r,w,s,1);
end

%% Visualize the vector field and apply independent transformations
if linear_registration==0
    flow = opticalFlow(-1*squeeze(disp_field(:,:,1)),-1*squeeze(disp_field(:,:,2)));
%     figure(), subplot(121), plot(flow, 'DecimationFactor', [8,16]), title('Overall Displacement Field');
%     subplot(122), imshow(abs_disp_field, []), colormap('jet'), title('Absolute Disp Field')

%% Warping only the frequency axis
    mod_disp_field = disp_field;
    mod_disp_field(:,:,1) = zeros(size(squeeze(mod_disp_field(:,:,1))));
    abs_warp = squeeze(mod_disp_field(:,:,1)).^2 + squeeze(mod_disp_field(:,:,2)).^2;

    registered_mag = imwarp(registered_mag,mod_disp_field);
    registered_phase = imwarp(registered_phase,mod_disp_field);

    flow = opticalFlow(-1*squeeze(mod_disp_field(:,:,1)),-1*squeeze(mod_disp_field(:,:,2)));
%     figure(), subplot(121), plot(flow, 'DecimationFactor',[8,16]), title('Freq Axis Warping field');
%     subplot(122), imshow(abs_warp, []), colormap('jet'), title('Absolute Warping field');

%     figure(), subplot(121), imshowpair(res_spect_neu_phase,res_spect_ang_phase), title('Unregistered');
%     subplot(122), imshowpair(registered_phase,res_spect_ang_phase), title('Freq Phase Registration');

%     figure(), subplot(121), imshowpair(res_spect_neu_mag,res_spect_ang_mag), title('Unregistered');
%     subplot(122), imshowpair(registered_mag,res_spect_ang_mag), title('Freq Spectrogram Registration');

%     figure(), subplot(121), imshow(spect_neu_mag,[]), colormap('jet'), title('Original Spectrogram');
%     subplot(122), imshow(registered_mag,[]), colormap('jet'), title('Freq Modified Spectrogram');

%     res_reg_mag = imresize(registered_mag,size(spect_neu_mag));
%     res_reg_phase = imresize(registered_phase,size(spect_neu_phase));
    recon_speech_freq_warped = get_speech(registered_mag,registered_phase,f,r,w,s,1);

%% Warping only the time axis
    mod_disp_field = disp_field;
    mod_disp_field(:,:,2) = zeros(size(squeeze(mod_disp_field(:,:,2))));
    abs_warp = squeeze(mod_disp_field(:,:,1)).^2 + squeeze(mod_disp_field(:,:,2)).^2;

    registered_mag = imwarp(registered_mag,mod_disp_field);
    registered_phase = imwarp(registered_phase,mod_disp_field);

    flow = opticalFlow(-1*squeeze(mod_disp_field(:,:,1)),-1*squeeze(mod_disp_field(:,:,2)));
%     figure(), subplot(121), plot(flow, 'DecimationFactor',[8,16]), title('Time Axis Warping field');
%     subplot(122), imshow(abs_warp, []), colormap('jet'), title('Absolute Warping field');

%     figure(), subplot(121), imshowpair(res_spect_neu_phase,res_spect_ang_phase), title('Unregistered');
%     subplot(122), imshowpair(registered_phase,res_spect_ang_phase), title('Time Phase Registration');

%     figure(), subplot(121), imshowpair(res_spect_neu_mag,res_spect_ang_mag), title('Unregistered');
%     subplot(122), imshowpair(registered_mag,res_spect_ang_mag), title('Time Spectrogram Registration');

%     figure(), subplot(121), imshow(spect_neu_mag,[]), colormap('jet'), title('Original Spectrogram');
%     subplot(122), imshow(registered_mag,[]), colormap('jet'), title('Time Modified Spectrogram');

%     res_reg_mag = imresize(registered_mag,size(spect_neu_mag));
%     res_reg_phase = imresize(registered_phase,size(spect_neu_phase));
    recon_speech_time_warped = get_speech(registered_mag,spect_ang_phase,f,r,w,s,1);
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