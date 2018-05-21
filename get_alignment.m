function [d1_recon, d2_recon] = get_alignment(d1,d2,sr)
    addpath(genpath('./dtw'));
    addpath(genpath('./mfcc'));

    d1 = -1 + 2*((d1 - min(d1)) / (max(d1) - min(d1)));
    d2 = -1 + 2*((d2 - min(d2)) / (max(d2) - min(d2)));

    [d1_mag, d1_phase] = get_spectrogram(d1,sr,512,0.025,0.015,0);
    [d2_mag, d2_phase] = get_spectrogram(d2,sr,512,0.025,0.015,0);

    [d1_feat,~,~] = mfcc(d1,sr,25,15,0.97,'hamming',[130, 4000],20,26,22);
    [d2_feat,~,~] = mfcc(d2,sr,25,15,0.97,'hamming',[130, 4000],20,26,22);

    M = simmx(d1_feat, d2_feat, 'cosine');
%     figure(), imshow(M, [])
    [p,q,~] = dp(1-M);
%     hold on; plot(q,p,'r'); hold off

    d1_mag_tilda = d1_mag(:,p);
    d1_phase_tilda = d1_phase(:,p);

    d2_mag_tilda = d2_mag(:,q);
    d2_phase_tilda = d2_phase(:,q);

%     figure(), subplot(121), imshow(d1_mag, []), title('D1 Original'), subplot(122), imshow(d1_mag_tilda, []), title('D1 Warped'), colormap('jet');
%     figure(), subplot(121), imshow(d2_mag, []), title('D2 Original'), subplot(122), imshow(d2_mag_tilda, []), title('D2 Warped'), colormap('jet');
    d1_recon = get_speech(d1_mag_tilda,d1_phase_tilda,sr,512,0.025,0.015,1);
    d2_recon = get_speech(d2_mag_tilda,d2_phase_tilda,sr,512,0.025,0.015,1);
%     sound(d1_recon, sr);
%     pause(2)
%     sound(d2_recon, sr);
end