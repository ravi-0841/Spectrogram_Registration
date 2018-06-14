function [vx,vy] = my_diffeomorphic_demons(F, M, opts)
    
    if nargin<3;                         opts                 = struct();       end
    if ~isfield(opts,'only_freq');       opts.only_freq       = 0;              end
    if ~isfield(opts,'sigma_fluid');     opts.sigma_fluid     = 0.7;            end
    if ~isfield(opts,'sigma_diff');      opts.sigma_diff      = 1.0;            end
    if ~isfield(opts,'sigma_x');         opts.sigma_x         = 1.0;            end
    if ~isfield(opts,'sigma_i');         opts.sigma_i         = 1.0;            end
    if ~isfield(opts,'compositive');     opts.compositive     = 1;              end
    if ~isfield(opts,'max_iter');        opts.max_iter        = 500;            end
    if ~isfield(opts,'step');            opts.step            = 1.0;            end
    if ~isfield(opts,'stop_criterion');  opts.stop_criterion  = 0.01;           end
    if ~isfield(opts,'pyramid_levels');  opts.pyramid_levels  = 2;              end
    
    vx = zeros(size(M));
    vy = zeros(size(M));
    
    M = mat2gray(M);
    F = mat2gray(F);
    
    for k = opts.pyramid_levels:-1:1
        scale_factor = 2^(-1*(k-1));
        
        M_tilda = imresize(M,scale_factor);
        F_tilda = imresize(F,scale_factor);
        
        vx = imresize(vx*scale_factor,scale_factor);
        vy = imresize(vy*scale_factor,scale_factor);
        
        [vx,vy] = register_images(F_tilda, M_tilda, vx, vy, opts);
        
        vx = imresize(vx/scale_factor, size(M));
        vy = imresize(vy/scale_factor, size(M));
    end
end