files = dir('./wav/*.wav');
I_cell = cell(length(files),1);

N = 1024;
wshift = 128;

W = hann(N);
S = hann(N);

for i = 1:length(files)
    [x,fs] = audioread(files(i).name);
    X = stft(x,N,wshift,W);
    X_mag = abs(X);
    X_ang = angle(X);
    
    disp(files(i).name);
    
    %% Deconvolution to get discrete objects
    PSF = fspecial('gaussian', [5 5], 5.0); % [5 5] works best
    I_trans = mat2gray(log(1+X_mag));
    I_trans = imdiffusefilt(I_trans); % Not sure about this step although it seems to be working
    I_decon = deconvlucy(I_trans, PSF, 10);
    I_cell{i,1} = I_decon;
    
    level = graythresh(I_decon);
    I_thresh = imbinarize(I_decon, level);
    ero_strel = [1 1; 1 1];
    I_eroded = imerode(I_thresh, ero_strel);
%     I_eroded = imdilate(I_eroded, [1 1]);

    subplot(131), imshow(I_trans,[]), subplot(132), imshow(I_decon, []), colormap(jet), ...
        subplot(133), imshow(I_eroded,[]);
    pause;
end