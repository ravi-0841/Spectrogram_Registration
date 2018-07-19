if ~ispc
    files       = dir('/home/ravi/Downloads/data_shape/*.wav');
    data_loc    = '/home/ravi/Downloads/data_shape/';
else
    files       = dir('C:\Users\Ravi Shankar\Downloads\data_shape\');
    data_loc    = 'C:\Users\Ravi Shankar\Downloads\data_shape\';
end
N               = 1024;
wshift          = 128;

W = hann(N);
S = hann(N);

%% 
ero_strel       = [1;1];
PSF             = fspecial('gaussian', [5 5], 5.0);
min_size        = 30;
num_coeffs      = 256;

shapes_angry    = [];
shapes_happy    = [];
shapes_nutrl    = [];

fc              = 1000;
y               = 0:N/2;
z               = 1 ./ sqrt(1 + (y./(fc*(N+2)/16000)).^10);

rand_idx        = randperm(floor(length(files)/3));

for iter       = 1:length(rand_idx)
    i          = rand_idx(iter);
    try
        [xa,~] = audioread([data_loc 'angry' num2str(i) '.wav']);
        [xh,~] = audioread([data_loc 'happy' num2str(i) '.wav']);
        [xn,~] = audioread([data_loc 'neutral' num2str(i) '.wav']);
    catch
        disp(['Missed i = ' num2str(i)]);
    end
    
    Xa = stft(xa,N,wshift,W);
    Xh = stft(xh,N,wshift,W);
    Xn = stft(xn,N,wshift,W);
    
    Xa_mag = abs(Xa);
    Xh_mag = abs(Xh);
    Xn_mag = abs(Xn);

    Xa_mag = Xa_mag.*(z' * ones(1,size(Xa_mag,2)));
    Xh_mag = Xh_mag.*(z' * ones(1,size(Xh_mag,2)));
    Xn_mag = Xn_mag.*(z' * ones(1,size(Xn_mag,2)));
    
    % Pre-processing to get discrete objects
    Ia_trans        = mat2gray(log(1+Xa_mag));
    Ia_trans        = imdiffusefilt(Ia_trans);
    Ia_decon        = deconvlucy(Ia_trans, PSF, 10);
    level           = graythresh(Ia_decon);
    Ia_thresh       = imbinarize(Ia_decon, level);
    Ia_eroded       = imopen(Ia_thresh, ero_strel);
    Ia_eroded       = mask_thickening(Ia_eroded);
    [~,angl]        = imgradient(Ia_eroded);
    ver_a           = (angl >= 170 | angl <= -170);
    ver_a           = circshift(ver_a, -1, 2);
    ver_a           = bwareaopen(ver_a, 3);
    Ia_eroded       = Ia_eroded - ver_a;
    Ia_eroded       = bwareaopen(Ia_eroded, min_size, 4);
    shapes_angry    = [shapes_angry get_fourier_descriptors(Ia_eroded,num_coeffs)];
    
    Ih_trans        = mat2gray(log(1+Xh_mag));
    Ih_trans        = imdiffusefilt(Ih_trans);
    Ih_decon        = deconvlucy(Ih_trans, PSF, 10);
    level           = graythresh(Ih_decon);
    Ih_thresh       = imbinarize(Ih_decon, level);
    Ih_eroded       = imopen(Ih_thresh, ero_strel);
    Ih_eroded       = mask_thickening(Ih_eroded);
    [~,angl]        = imgradient(Ih_eroded);
    ver_h           = (angl >= 170 | angl <= -170);
    ver_h           = bwareaopen(ver_h, 3);
    ver_h           = circshift(ver_h, -1, 2);
    Ih_eroded       = Ih_eroded - ver_h;
    Ih_eroded       = bwareaopen(Ih_eroded, min_size, 4);
    shapes_happy    = [shapes_happy get_fourier_descriptors(Ih_eroded,num_coeffs)];
    
    In_trans        = mat2gray(log(1+Xn_mag));
    In_trans        = imdiffusefilt(In_trans);
    In_decon        = deconvlucy(In_trans, PSF, 10);
    level           = graythresh(In_decon);
    In_thresh       = imbinarize(In_decon, level);
    In_eroded       = imopen(In_thresh, ero_strel);
    In_eroded       = mask_thickening(In_eroded);
    [~,angl]        = imgradient(In_eroded);
    ver_n           = (angl >= 170 | angl <= -170);
    ver_n           = bwareaopen(ver_n, 3);
    ver_n           = circshift(ver_n, -1, 2);
    In_eroded       = In_eroded - ver_n;
    In_eroded       = bwareaopen(In_eroded, min_size, 4);
    shapes_nutrl    = [shapes_nutrl get_fourier_descriptors(In_eroded,num_coeffs)];
    
%     figure(1), imshow(Ia_eroded,[]), title('Angry');
%     figure(2), imshow(Ih_eroded,[]), title('Happy');
%     figure(3), imshow(In_eroded,[]), title('Neutral');
%     pause;
end

clear files N wshift W S ero_strel PSF min_size num_coeffs fc y z Ia_trans Ih_trans In_trans ...
    xa xh xn Xa Xn Xh Xa_mag Xh_mag Xn_mag Ia_decon Ih_decon In_decon Ia_thresh Ih_thresh In_thresh ...
    ver_a ver_h ver_n level Ia_eroded Ih_eroded In_eroded angl i iter rand_idx

%% KSVD shapes
load('./data/dict_angry.mat');
load('./data/dict_happy.mat');
load('./data/dict_neutral.mat');

for i = 1:50
    z1 = dict_ang(:,i); 
    z2 = dict_hap(:,i); 
    z3 = dict_nut(:,i);
    figure(1),...
        subplot(131),plotEllipticFourierDescriptor(z1(1:256),z1(257:512),z1(513:768),z1(769:1024),100,1),...
        title('Angry'),subplot(132),...
        plotEllipticFourierDescriptor(z2(1:256),z2(257:512),z2(513:768),z2(769:1024),100,1),...
        title('Happy'),subplot(133),...
        plotEllipticFourierDescriptor(z3(1:256),z3(257:512),z3(513:768),z3(769:1024),100,1),title('Neutral');
    pause;
end

%% PCA
if ~ispc
    load('/home/ravi/Downloads/256_fourier_shape_data.mat');
else
    load('C:\Users\Ravi Shankar\Downloads\256_fourier_shape_data.mat');
end
[ang_coeff, ang_score, ang_latent] = pca(shapes_angry');
[hap_coeff, hap_score, hap_latent] = pca(shapes_happy');
[neu_coeff, neu_score, neu_latent] = pca(shapes_nutrl');

for i = 1:50
    z1 = ang_coeff(:,i); 
    z2 = hap_coeff(:,i); 
    z3 = neu_coeff(:,i);
    figure(1),...
        subplot(131),plotEllipticFourierDescriptor(z1(1:256),z1(257:512),z1(513:768),z1(769:1024),100,1),...
        title('Angry'),subplot(132),...
        plotEllipticFourierDescriptor(z2(1:256),z2(257:512),z2(513:768),z2(769:1024),100,1),...
        title('Happy'),subplot(133),...
        plotEllipticFourierDescriptor(z3(1:256),z3(257:512),z3(513:768),z3(769:1024),100,1),title('Neutral');
    pause;
end

complete_data = [shapes_angry shapes_happy shapes_nutrl];
[full_coeff, full_score, full_latent] = pca(complete_data');

figure(2);
scatter3(full_score(1:length(shapes_angry),1),...
    full_score(1:length(shapes_angry),2),full_score(1:length(shapes_angry),3), 'r');
hold on, scatter3(full_score(length(shapes_angry)+1:length(shapes_angry)+length(shapes_happy),1),...
    full_score(length(shapes_angry)+1:length(shapes_angry)+length(shapes_happy),2),...
    full_score(length(shapes_angry)+1:length(shapes_angry)+length(shapes_happy),3), 'g');
hold on, scatter3(full_score(length(shapes_angry)+length(shapes_happy)+1:end,1),...
    full_score(length(shapes_angry)+length(shapes_happy)+1:end,2),...
    full_score(length(shapes_angry)+length(shapes_happy)+1:end,3), 'b');

figure(3);
plot(full_score(1:length(shapes_angry),1),...
    full_score(1:length(shapes_angry),2), 'ro');
hold on, plot(full_score(length(shapes_angry)+1:length(shapes_angry)+length(shapes_happy),1),...
    full_score(length(shapes_angry)+1:length(shapes_angry)+length(shapes_happy),2), 'go');
hold on, plot(full_score(length(shapes_angry)+length(shapes_happy)+1:end,1),...
    full_score(length(shapes_angry)+length(shapes_happy)+1:end,2), 'bo');

% q = randperm(min([length(shapes_angry),length(shapes_happy),length(shapes_nutrl)]));
% figure(4)
% scatter3(ang_score(q(1:100),1),ang_score(q(1:100),2),ang_score(q(1:100),3), 'r');
% hold on, scatter3(hap_score(q(1:100),1),hap_score(q(1:100),2),hap_score(q(1:100),3), 'g');
% hold on, scatter3(neu_score(q(1:100),1),neu_score(q(1:100),2),neu_score(q(1:100),3), 'b');