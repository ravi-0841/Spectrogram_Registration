function [d1_recon, d2_recon] = get_alignment(d1,d2,f,w,o,r,topo)
    addpath(genpath('./dtw'));
    addpath(genpath('./mfcc'));
    
    monotonicity = {'unconstrained', 'slopeHalf', 'slopeThird'};
    s = w - o;

    d1 = -1 + 2*((d1 - min(d1)) / (max(d1) - min(d1)));
    d2 = -1 + 2*((d2 - min(d2)) / (max(d2) - min(d2)));

    d1_spect = spectrogram(d1,w*f,int64(o*f),r);
    d2_spect = spectrogram(d2,w*f,int64(o*f),r);
    
    d1_mag = abs(d1_spect);
    d1_phase = angle(d1_spect);
    
    d2_mag = abs(d2_spect);
    d2_phase = angle(d2_spect);

    [d1_feat,~,~] = mfcc(d1,f,w*1000,s*1000,0.97,'hamming',[130, 4000],20,26,22);
    [d2_feat,~,~] = mfcc(d2,f,w*1000,s*1000,0.97,'hamming',[130, 4000],20,26,22);
    
    M = pair_sim(d1_feat,d2_feat,'cosine');
    [p,q,~] = dtw(1-M,monotonicity{topo});
%     figure(), subplot(121), imshow(D, []), title(monotonicity{topo})
%     hold on; plot(q,p,'r'); hold off
    
    %%

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

%     figure(), subplot(121), imshow(d1_mag, []), title('D1 Original'), subplot(122), imshow(d1_mag_tilda, []), title('D1 Warped'), colormap('jet');
%     figure(), subplot(121), imshow(d2_mag, []), title('D2 Original'), subplot(122), imshow(d2_mag_tilda, []), title('D2 Warped'), colormap('jet');
    d1_recon = get_speech(d1_mag_tilda,d1_phase_tilda,f,r,w,s,1);
    d2_recon = get_speech(d2_mag_tilda,d2_phase_tilda,f,r,w,s,1);
%     sound(d1_recon, f);
%     pause(2)
%     sound(d2_recon, f);
end