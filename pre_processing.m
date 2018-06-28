files = dir('./wav/*.wav');
I_cell = cell(length(files),1);

N = 1024;
wshift = 128;
epsilon = 0;

W = hann(N);
S = hann(N);

resize_scale = 1;

% for i = 1:length(files)
%     [x,fs] = audioread(files(i).name);
%     X = stft(x,N,wshift,W);
%     X_mag = abs(X);
%     X_ang = angle(X);
%     
%     disp(files(i).name);
%     
%     % Deconvolution to get discrete objects
%     PSF = fspecial('gaussian', resize_scale*[5 5], resize_scale*5.0); % [5 5] and 5.0 works best for N = 1024
%     I_trans = imresize(mat2gray(log(1+X_mag)), resize_scale);
%     I_sharp = imsharpen(I_trans);
%     I_trans = imdiffusefilt(I_trans); % Not sure about this step
%     I_decon = deconvlucy(I_trans, PSF, 10);
%     
%     level = graythresh(I_decon);
%     I_thresh = imbinarize(I_decon, level);
%     
% %     level = adaptthresh(I_decon);
% %     I_thresh = imbinarize(I_decon, level);
% 
% %     idx = kmeans(I_decon(:), 2);
% %     I_thresh = zeros(size(I_decon));
% %     intensity1 = mean(mean(I_decon(idx==1)));
% %     intensity2 = mean(mean(I_decon(idx==2)));
% %     if intensity1 > intensity2
% %         I_thresh(idx==1) = 1;
% %     else
% %         I_thresh(idx==2) = 1;
% %     end
%     
% %     level = otsuthresh(imhist(I_decon));
% %     I_thresh = imbinarize(I_decon, level);
%     
%     ero_strel = [1;1]; % So far this structuring element has done OK
%     I_eroded = imopen(I_thresh, ero_strel);
%     I_cell{i,1} = I_eroded;
%     I_eroded = mask_thickening(I_eroded);
%     I_eroded = bwareaopen(I_eroded,30);
%     
%     connected_objs = bwconncomp(I_eroded, 4);
% 
%     subplot(131), imshow(I_trans,[]), subplot(132), imshow(I_decon, []), colormap(jet), ...
%         subplot(133), imshow(I_eroded,[]);
%     pause;
%     
% end

%% 
ero_strel = [1;1]; % So far this structuring element has done OK
PSF = fspecial('gaussian', resize_scale*[5 5], resize_scale*5.0); % [5 5] and 5.0 works best for N = 1024

for i = 1:6
    [xa,fs] = audioread(['angry' num2str(i) '.wav']);
    [xh,fs] = audioread(['happy' num2str(i) '.wav']);
    [xn,fs] = audioread(['neutral' num2str(i) '.wav']);
    Xa = stft(xa,N,wshift,W);
    Xh = stft(xh,N,wshift,W);
    Xn = stft(xn,N,wshift,W);
    
    Xa_mag = abs(Xa);
    Xa_ang = angle(Xa);
    
    Xh_mag = abs(Xh);
    Xh_ang = angle(Xh);
    
    Xn_mag = abs(Xn);
    Xn_ang = angle(Xn);
    
    % Deconvolution to get discrete objects
    Ia_trans = imresize(mat2gray(log(1+Xa_mag)), resize_scale);
    Ia_trans = imdiffusefilt(Ia_trans);
    Ia_decon = deconvlucy(Ia_trans, PSF, 10);
    level = graythresh(Ia_decon);
    Ia_thresh = imbinarize(Ia_decon, level-epsilon);
    Ia_eroded = imopen(Ia_thresh, ero_strel);
    Ia_eroded = mask_thickening(Ia_eroded);
    Ia_eroded = bwareaopen(Ia_eroded, 30);
    rgpa = regionprops(bwconncomp(Ia_eroded, 4), 'centroid');
    centroid_a = cat(1, rgpa.Centroid);
    
    Ih_trans = imresize(mat2gray(log(1+Xh_mag)), resize_scale);
    Ih_trans = imdiffusefilt(Ih_trans);
    Ih_decon = deconvlucy(Ih_trans, PSF, 10);
    level = graythresh(Ih_decon);
    Ih_thresh = imbinarize(Ih_decon, level-epsilon);
    Ih_eroded = imopen(Ih_thresh, ero_strel);
    Ih_eroded = mask_thickening(Ih_eroded);
    Ih_eroded = bwareaopen(Ih_eroded, 30);
    rgph = regionprops(bwconncomp(Ih_eroded, 4), 'centroid');
    centroid_h = cat(1, rgph.Centroid);
    
    In_trans = imresize(mat2gray(log(1+Xn_mag)), resize_scale);
    In_trans = imdiffusefilt(In_trans);
    In_decon = deconvlucy(In_trans, PSF, 10);
    level = graythresh(In_decon);
    In_thresh = imbinarize(In_decon, level-epsilon);
    In_eroded = imopen(In_thresh, ero_strel);
    In_eroded = mask_thickening(In_eroded);
    In_eroded = bwareaopen(In_eroded, 30);
    rgpn = regionprops(bwconncomp(In_eroded, 4),'centroid');
    centroid_n = cat(1, rgpn.Centroid);
    
    subplot(131), imshow(Ia_eroded.*Ia_trans,[]), hold on, plot(centroid_a(:,1), centroid_a(:,2), 'g.'), hold off, ...
        title('Angry'), subplot(132), imshow(Ih_eroded.*Ih_trans, []), hold on, plot(centroid_h(:,1), centroid_h(:,2), 'g.'), hold off, ...
        title('Happy'), subplot(133), imshow(In_eroded.*In_trans,[]), hold on, plot(centroid_n(:,1), centroid_n(:,2), 'g.'), hold off, title('Neutral');
    pause;
end