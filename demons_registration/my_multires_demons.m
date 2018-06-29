function disp_field = my_multires_demons(F, M, N, opts)
    
    if nargin<4
        opts = struct();
    end
    
    if nargin<3
        N = 512;
    end
    
    if ~isfield(opts,'only_freq');       opts.only_freq       = 0;              end
    if ~isfield(opts,'sigma_fluid');     opts.sigma_fluid     = 0.7;            end
    if ~isfield(opts,'sigma_diff');      opts.sigma_diff      = 0.7;            end
    if ~isfield(opts,'sigma_x');         opts.sigma_x         = 1.0;            end
    if ~isfield(opts,'sigma_i');         opts.sigma_i         = 1.0;            end
    if ~isfield(opts,'lambda');          opts.lambda          = 0.5;            end
    if ~isfield(opts,'window_size');     opts.window_size     = 50;             end
    if ~isfield(opts,'stride');          opts.stride          = 50;             end
    if ~isfield(opts,'compositive');     opts.compositive     = 1;              end
    if ~isfield(opts,'max_iter');        opts.max_iter        = 500;            end
    if ~isfield(opts,'max_epochs');      opts.max_epochs      = 1;              end
    if ~isfield(opts,'step');            opts.step            = 1.0;            end
    if ~isfield(opts,'stop_criterion');  opts.stop_criterion  = 0.01;           end
    if ~isfield(opts,'pyramid_levels');  opts.pyramid_levels  = 2;              end
    if ~isfield(opts,'diffeomorphism');  opts.diffeomorphism  = 0;              end
    if ~isfield(opts,'plot');            opts.plot            = 0;              end
    
%     I_cell = cell(opts.pyramid_levels, opts.max_iter);
    
    vx = zeros(size(M));
    vy = zeros(size(M));
    
    M = mat2gray(M);
    F = mat2gray(F);
    
    if opts.plot
        figure(1)
        subplot(121), title('Difference');
        subplot(122), title('SSD');
    end
    
    for k = opts.pyramid_levels:-1:1
        scale_factor = 2^(-1*(k-1));
        
        M_tilda = imresize(M,scale_factor);
        F_tilda = imresize(F,scale_factor);
        
        vx = imresize(vx*scale_factor,scale_factor);
        vy = imresize(vy*scale_factor,scale_factor);
        
        if opts.diffeomorphism
            disp_field = my_diffeomorphism(F_tilda, M_tilda, N, vx, vy, opts);
        else
            disp_field = window_demons_hoz(F_tilda, M_tilda, vx, vy, opts); % change to my_demons for original demons
        end
        
%         for image_num = 1:opts.max_iter
%             I_cell{k,image_num} = imresize(cell_img{image_num,1}, [200*(3-k),150*(3-k)]);
%         end
        
        vx = squeeze(disp_field(:,:,1));
        vy = squeeze(disp_field(:,:,2));
        
        vx = imresize(vx/scale_factor, size(M));
        vy = imresize(vy/scale_factor, size(M));
        
        if opts.plot;   subplot(122), hold on; end
    end
end