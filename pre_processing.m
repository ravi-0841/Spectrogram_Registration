files = dir('./wav/*.wav');
I_cell = cell(length(files),1);

N = 1024;
wshift = 128;

W = hann(N);
S = hann(N);

resize_scale = 1;

for i = 1:length(files)
    [x,fs] = audioread(files(i).name);
    X = stft(x,N,wshift,W);
    X_mag = abs(X);
    X_ang = angle(X);
    
    disp(files(i).name);
    
    %% Deconvolution to get discrete objects
    PSF = fspecial('gaussian', resize_scale*[5 5], resize_scale*5.0); % [5 5] and 5.0 works best for N = 1024
    I_trans = imresize(mat2gray(log(1+X_mag)), resize_scale);
    I_sharp = imsharpen(I_trans);
    I_trans = imdiffusefilt(I_trans); % Not sure about this step although it seems to be working
    I_decon = deconvlucy(I_trans, PSF, 10);
    I_cell{i,1} = I_decon;
    
    level = graythresh(I_decon);
    I_thresh = imbinarize(I_decon, level);
    ero_strel = [1;1]; % So far this structuring element has done OK
    I_eroded = imopen(I_thresh, ero_strel);

    subplot(131), imshow(I_trans,[]), subplot(132), imshow(I_decon, []), colormap(jet), ...
        subplot(133), imshow(I_eroded,[]);
    pause;
end
