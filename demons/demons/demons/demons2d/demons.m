% Demons Registration
%
% Simple matlab code for 2D image registration using the diffeomorphic log-demons algorithm 
% Code is provided in order to help the understanding of the Demons algorithm
%
% Original Algorithm in:
% [1] - Symmetric Log-Domain Diffeomorphic Registration: A Demons-Based Approach
%       Vercauteren, Pennec, Perchant, Ayache -- MICCAI 2008, 754-761
% [2] - Diffeomorphic demons: Efficient non-parametric image registration,
%       Vercauteren, Pennec, Perchant, Ayache -- NeuroImage 2009, (45)1:61-72
%
% For a more recent work/survey, exploiting global shape characteristics 
% (instead of the conventional local gradient-based approaches), consider citing 
%
% [1] - Spectral Log-Demons: Diffeomorphic Image Registration with Very Large Deformations
%       Lombaert, Grady, Pennec, Ayache, Cheriet -- IJCV 2014, (107)3:254-271


function demons(F,M)

figure(1); clf; colormap gray; set(gcf,'renderer','painter')

%% Parameters
niter           = 500;
sigma_fluid     = 2.0; % regularize update      field
sigma_diffusion = 2.0; % regularize deformation field
sigma_i         = 1.0; % weight on similarity term
sigma_x         = 1.0; % weight on spatial uncertainties (maximal step)
nlevel          = 3;   % multiresolution
do_display      = 1;   % display iterations

%% Load fixed image
% F = double(imread('data/lenag2.png'));          % fixed   image
%F = double(imread('statue-rio.png'));           % fixed   image
F = imresize(F,0.5);
F = 255*(F-min(F(:)))/range(F(:));              % normalize intensities

%% Load moving image
% M = double(imread('data/lenag1.png'));          % moving  image
%M = double(imread('statue-rio-deformed.png'));  % moving  image
M = imresize(M,0.5);
M = 255*(M-min(M(:)))/range(M(:));              % normalize intensities

% Translate
%shift = 3; tmp = zeros(size(M)); tmp((1+shift):end,:) = M(1:(end-shift),:); M = tmp;

%% Create random moving image
%[F,sx0,sy0]  = randomdeform(F,50,5);
%figure(5); showvector(sx0,sy0,4,3,lim); drawnow;

if nlevel == 1
    %% Register
    disp(['Register...']);
    opt = struct('niter',niter, 'sigma_fluid',sigma_fluid, 'sigma_diffusion',sigma_diffusion, 'sigma_i',sigma_i, 'sigma_x',sigma_x, 'do_display',do_display, 'do_plotenergy',1);
    [Mp,sx,sy,vx,vy] = register(F,M,opt);

else 
    %% Multiresolution
    vx = zeros(size(M)); % deformation field
    vy = zeros(size(M));
    for k=nlevel:-1:1
        disp(['Register level: ' num2str(k) '...']);

        % downsample
        scale = 2^-(k-1);
        Fl = imresize(F,scale);
        Ml = imresize(M,scale);
        vxl = imresize(vx*scale,scale);
        vyl = imresize(vy*scale,scale);

        % register
        opt = struct('niter',niter,...
            'sigma_fluid',sigma_fluid,...
            'sigma_diffusion',sigma_diffusion,...
            'sigma_i',sigma_i,...
            'sigma_x',sigma_x,...
            'vx',vxl, 'vy',vyl,...
            'do_display',do_display, 'do_plotenergy',1);
        [Mp,sxl,syl,vxl,vyl] = register(Fl,Ml,opt);

        % upsample
        vx = imresize(vxl/scale,size(M));
        vy = imresize(vyl/scale,size(M));
    end
    [sx,sy] = expfield(vx,vy);
    
end
end