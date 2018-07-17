files           = dir('./wav/*.wav');
N               = 1024;
wshift          = 128;

W = hann(N);
S = hann(N);

%% 
ero_strel       = [1;1];
PSF             = fspecial('gaussian', resize_scale*[5 5], resize_scale*5.0);
min_size        = 30;
num_coeffs      = 64;

shapes_angry    = [];
shapes_happy    = [];
shapes_nutrl    = [];

fc              = 1000;
y               = 0:N/2;
z               = 1 ./ sqrt(1 + (y./(fc*(N+2)/16000)).^10);

for i = 1:length(files)/3
    [xa,~] = audioread(['angry' num2str(i) '.wav']);
    [xh,~] = audioread(['happy' num2str(i) '.wav']);
    [xn,~] = audioread(['neutral' num2str(i) '.wav']);
    
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
    In_thresh       = imbinarize(In_decon, level-epsilon);
    In_eroded       = imopen(In_thresh, ero_strel);
    In_eroded       = mask_thickening(In_eroded);
    [~,angl]        = imgradient(In_eroded);
    ver_n           = (angl >= 170 | angl <= -170);
    ver_n           = bwareaopen(ver_n, 3);
    ver_n           = circshift(ver_n, -1, 2);
    In_eroded       = In_eroded - ver_n;
    In_eroded       = bwareaopen(In_eroded, min_size, 4);
    shapes_nutrl    = [shapes_nutrl get_fourier_descriptors(In_eroded,num_coeffs)];
end

clear files N wshift W S ero_strel PSF min_size num_coeffs fc y z Ia_trans Ih_trans In_trans ...
    xa xh xn Xa Xn Xh Xa_mag Xh_mag Xn_mag Ia_decon Ih_decon In_decon Ia_thresh Ih_thresh In_thresh ...
    ver_a ver_h ver_n level Ia_eroded Ih_eroded In_eroded angl