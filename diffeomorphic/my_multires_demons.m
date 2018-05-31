function disp_field = my_multires_demons(F, M, opts)
    
    if nargin<3;                         opts                 = struct();       end
    if ~isfield(opts,'sigma_fluid');     opts.sigma_fluid     = 0.7;            end
    if ~isfield(opts,'sigma_diff');      opts.sigma_diff      = 1.0;            end
    if ~isfield(opts,'sigma_x');         opts.sigma_x         = 1.0;            end
    if ~isfield(opts,'sigma_i');         opts.sigma_i         = 1.0;            end
    if ~isfield(opts,'compositive');     opts.compositive     = 1;              end
    if ~isfield(opts,'max_iter');        opts.max_iter        = 500;            end
    if ~isfield(opts,'step');            opts.step            = 1.0;            end
    if ~isfield(opts,'stop_criterion');  opts.stop_criterion  = 0.01;           end
    if ~isfield(opts,'pyramid_levels');  opts.pyramid_levels  = 1;              end
    if ~isfield(opts,'diffeomorphism');  opts.diffeomorphism  = 0;              end
    
    vx = zeros(size(M));
    vy = zeros(size(M));
    
    M = mat2gray(M);
    F = mat2gray(F);
    
    figure(1)
    subplot(121), title('Difference');
    subplot(122), title('SSD');
    
    for k = opts.pyramid_levels:-1:1
        scale_factor = 2^(-1*(k-1));
        
        M_tilda = imresize(M,scale_factor);
        F_tilda = imresize(F,scale_factor);
        
        vx = imresize(vx*scale_factor,scale_factor);
        vy = imresize(vy*scale_factor,scale_factor);
        
        if opts.diffeomorphism
            disp_field = my_diffeomorphism(F_tilda, M_tilda, vx, vy, opts);
        else
            disp_field = my_demons(F_tilda, M_tilda, vx, vy, opts);
        end
        
        vx = squeeze(disp_field(:,:,1));
        vy = squeeze(disp_field(:,:,2));
        
        vx = imresize(vx/scale_factor, size(M));
        vy = imresize(vy/scale_factor, size(M));
        subplot(122), hold on;
    end
end