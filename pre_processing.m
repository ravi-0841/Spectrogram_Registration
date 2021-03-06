files = dir('./wav/*.wav');
I_cell = cell(length(files),1);

N = 1024;
wshift = 128;
epsilon = 0;

W = hann(N);
S = hann(N);

resize_scale = 1;

r = 512;
w = 0.025;
s = 0.010;
top = 3;

ero_strel = [1;1];
PSF = fspecial('gaussian', resize_scale*[5 5], resize_scale*5.0);
min_size = 30;

for i = 1:floor(length(files)/4)
    try
        [xa,fs] = audioread(['angry' num2str(i) '.wav']);
        [xh,fs] = audioread(['happy' num2str(i) '.wav']);
        [xn,fs] = audioread(['neutral' num2str(i) '.wav']);
    catch
    end
    
%     [xn, xa] = get_alignment(xn,xa,fs,w,w-s,r,top);
%     [xa, xh] = get_alignment(xa,xh,fs,w,w-s,r,top);
    
    Xa = stft(xa,N,wshift,W);
    Xh = stft(xh,N,wshift,W);
    Xn = stft(xn,N,wshift,W);
    
    Xa_mag = abs(Xa);
    Xa_ang = angle(Xa);
    
    Xh_mag = abs(Xh);
    Xh_ang = angle(Xh);
    
    Xn_mag = abs(Xn);
    Xn_ang = angle(Xn);
    
%     Xa_mag = imresize(Xa_mag, size(Xn_mag));
%     Xh_mag = imresize(Xh_mag, size(Xn_mag));
    
    fc = 2000;
    y = 0:N/2;
    z = 1 ./ sqrt(1 + (y./(fc*(N+2)/fs)).^10);
    Xa_mag = Xa_mag.*(z' * ones(1,size(Xa_mag,2)));
    Xh_mag = Xh_mag.*(z' * ones(1,size(Xh_mag,2)));
    Xn_mag = Xn_mag.*(z' * ones(1,size(Xn_mag,2)));
    
    % Deconvolution to get discrete objects
    Ia_trans = imresize(mat2gray(log(1+Xa_mag)), resize_scale);
    Ia_trans = imdiffusefilt(Ia_trans);
    Ia_decon = deconvlucy(Ia_trans, PSF, 10);
    level = graythresh(Ia_decon);
    Ia_thresh = imbinarize(Ia_decon, level-epsilon);
    Ia_eroded = imopen(Ia_thresh, ero_strel);
    Ia_eroded = mask_thickening(Ia_eroded);
    for line_removal = 1:3
        [~,angl]  = imgradient(Ia_eroded);
        ver_a     = (angl >= 170 | angl <= -170);
        ver_a     = circshift(ver_a, -1, 2);
        ver_a     = bwareaopen(ver_a,3);
        Ia_eroded = Ia_eroded - ver_a;
    end
    Ia_eroded = bwareaopen(Ia_eroded, min_size, 4);
    rgpa = regionprops(bwconncomp(Ia_eroded, 4), 'centroid');
    centroid_a = cat(1, rgpa.Centroid);
    
    Ih_trans = imresize(mat2gray(log(1+Xh_mag)), resize_scale);
    Ih_trans = imdiffusefilt(Ih_trans);
    Ih_decon = deconvlucy(Ih_trans, PSF, 10);
    level = graythresh(Ih_decon);
    Ih_thresh = imbinarize(Ih_decon, level-epsilon);
    Ih_eroded = imopen(Ih_thresh, ero_strel);
    Ih_eroded = mask_thickening(Ih_eroded);
    for line_removal = 1:3
        [~,angl]  = imgradient(Ih_eroded);
        ver_h     = (angl >= 170 | angl <= -170);
        ver_h     = bwareaopen(ver_h, 3);
        ver_h     = circshift(ver_h, -1, 2);
        Ih_eroded = Ih_eroded - ver_h;
    end
    Ih_eroded = bwareaopen(Ih_eroded, min_size, 4);
    rgph = regionprops(bwconncomp(Ih_eroded, 4), 'centroid');
    centroid_h = cat(1, rgph.Centroid);
    
    In_trans = imresize(mat2gray(log(1+Xn_mag)), resize_scale);
    In_trans = imdiffusefilt(In_trans);
    In_decon = deconvlucy(In_trans, PSF, 10);
    level = graythresh(In_decon);
    In_thresh = imbinarize(In_decon, level-epsilon);
    In_eroded = imopen(In_thresh, ero_strel);
    In_eroded = mask_thickening(In_eroded);
    for line_removal = 1:3
        [~,angl]    = imgradient(In_eroded);
        ver_n       = (angl >= 170 | angl <= -170);
        ver_n       = bwareaopen(ver_n, 3);
        ver_n       = circshift(ver_n, -1, 2);
        In_eroded   = In_eroded - ver_n;
    end
    In_eroded = bwareaopen(In_eroded, min_size, 4);
    rgpn = regionprops(bwconncomp(In_eroded, 4),'centroid');
    centroid_n = cat(1, rgpn.Centroid);
    
    %% Manual Fixing of the binary images     
    bwpaint1(Ia_eroded,Ia_decon);
    pause;
    manual_paint = load('manual_paint.mat');
    Ia_eroded = manual_paint.manual_paint;
    close all
    
    bwpaint1(Ih_eroded,Ih_decon);
    pause;
    manual_paint = load('manual_paint.mat');
    Ih_eroded = manual_paint.manual_paint;
    close all
    
    bwpaint1(In_eroded,In_decon);
    pause;
    manual_paint = load('manual_paint.mat');
    In_eroded = manual_paint.manual_paint;
    close all
    
    %% Skeletonization
    Ia_skeleton = bwskel(Ia_eroded);
    Ih_skeleton = bwskel(Ih_eroded);
    In_skeleton = bwskel(In_eroded);
    
%     [~,angl] = imgradient(Ia_skeleton);
%     hoz_a = (angl <= 100 | angl >= 80);
%     hoz_a = bwareaopen(hoz_a, 3);
    figure(1)
    subplot(131), imshow(Ia_skeleton), subplot(132), imshow(Ih_skeleton), subplot(133), imshow(In_skeleton);
%     imshowpair(Ia_skeleton,hoz_a);
%     disp('random')
    
    %% Scaling Component-wise
%     va = vertical_segmentation(Ia_eroded);
%     vh = vertical_segmentation(Ih_eroded);
%     vn = vertical_segmentation(In_eroded);
%     
%     bwa = bwconncomp(va);
%     bwh = bwconncomp(vh);
%     bwn = bwconncomp(vn);
%     
%     assert((bwa.NumObjects==bwh.NumObjects) && (bwa.NumObjects==bwn.NumObjects),...
%                                             'The utterances may not be same');
%     
%     new_im_a = zeros(size(Ia_eroded));
%     new_im_h = zeros(size(Ih_eroded));
%     new_im_n = zeros(size(In_eroded));
%     
%     for obj = 1:bwa.NumObjects
%         tempa = zeros(bwa.ImageSize);
%         temph = zeros(bwh.ImageSize);
%         tempn = zeros(bwn.ImageSize);
%         
%         tempa(bwa.PixelIdxList{obj}) = 1;
%         temph(bwh.PixelIdxList{obj}) = 1;
%         tempn(bwn.PixelIdxList{obj}) = 1;
%         
%         im_a_seg = Ia_eroded.*tempa;
%         im_h_seg = Ih_eroded.*temph;
%         im_n_seg = In_eroded.*tempn;
%         
%         [ra,~] = find(im_a_seg==1);
%         [rh,~] = find(im_h_seg==1);
%         [rn,~] = find(im_n_seg==1);
%         
%         ref = max([max(ra), max(rh), max(rn)]);
%         
%         im_a_seg = imresize(log(1+Xa_mag).*Ia_eroded, [513*ref/max(ra), size(im_a_seg,2)]);
%         im_h_seg = imresize(log(1+Xh_mag).*Ih_eroded, [513*ref/max(rh), size(im_h_seg,2)]);
%         im_n_seg = imresize(log(1+Xn_mag).*In_eroded, [513*ref/max(rn), size(im_n_seg,2)]);
%         
% %         im_a_seg(im_a_seg~=0) = 1;
% %         im_h_seg(im_h_seg~=0) = 1;
% %         im_n_seg(im_n_seg~=0) = 1;
%         
%         im_a_seg = im_a_seg(1:513,:);
%         im_h_seg = im_h_seg(1:513,:);
%         im_n_seg = im_n_seg(1:513,:);
%         
%         new_im_a(tempa==1) = im_a_seg(tempa==1);
%         new_im_h(temph==1) = im_h_seg(temph==1);
%         new_im_n(tempn==1) = im_n_seg(tempn==1);
%         
%     end

%     figure(2);
%     subplot(131), imshow(Ia_eroded,[]), hold on, plot(centroid_a(:,1), centroid_a(:,2), 'g.'), hold off, ...
%         title('Angry'), subplot(132), imshow(Ih_eroded,[]), hold on, ...
%         plot(centroid_h(:,1), centroid_h(:,2), 'g.'), hold off, ...
%         title('Happy'), subplot(133), imshow(In_eroded,[]), hold on, ...
%         plot(centroid_n(:,1), centroid_n(:,2), 'g.'), hold off, title('Neutral');
%     figure(3);
%     subplot(131), imshow(Ia_decon, []), subplot(132), imshow(Ih_decon, []),...
%       subplot(133), imshow(In_decon, []), colormap(jet);
    pause;
end






%%
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