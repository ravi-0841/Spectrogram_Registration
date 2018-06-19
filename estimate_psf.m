files = dir('./wav/*.wav');
I_cell = cell(length(files),1);

N = 512;
wshift = 128;

W = hann(N);
S = hann(N);

for i = 1:length(files)
    [x,fs] = audioread(files(i).name);
    X = stft(x,N,wshift,W);
    X_mag = abs(X);
    X_ang = angle(X);
    
    disp(files(i).name);
    
    %% Analysing Gradients
    I_trans = mat2gray(log(1+X_mag));
    I_gauss1 = imgaussfilt(I_trans, 0.8);
    I_gauss2 = imgaussfilt(I_trans, 1);
%     I_trans = imdiffusefilt(I_trans); % Not sure about this step although it seems to be working
    edge_detect = abs(I_gauss2 - I_gauss1);
    [Gx, Gy] = imgradientxy(I_trans);
    I_cell{i,1} = (Gx.^2 + Gy.^2);
    
    subplot(131), imshow(I_trans,[]), colormap(jet), subplot(132), imshow(Gx, []), ...
        subplot(133), imshow(edge_detect,[]);
    pause;
end
