function disp_field = my_diffeomorphism(F, M, N, vx, vy, opts)
    
    if nargin<3;                        opts                 = struct();     end
    if ~isfield(opts,'only_freq');      opts.only_freq       = 0;            end
    if ~isfield(opts,'alpha');          opts.alpha           = 0.4;          end
    if ~isfield(opts,'sigma_fluid');    opts.sigma_fluid     = 1.0;          end
    if ~isfield(opts,'sigma_diff');     opts.sigma_diff      = 1.0;          end
    if ~isfield(opts,'step');           opts.step            = 1.0;          end
    if ~isfield(opts,'epsilon');        opts.epsilon         = 10;           end
    if ~isfield(opts,'compositive');    opts.compositive     = 0;            end
    if ~isfield(opts,'max_iter');       opts.max_iter        = 1000;         end
    if ~isfield(opts,'plot');           opts.plot            = 0;            end
    

    disp(['Initial SSD ' num2str(sum(sum((F - M).^2)))]);
    disp(['Initial MI ' num2str(mutual_info(F, M))]);
    
    if isempty(vx) || isempty(vy)
        vec_field_x = zeros(size(M));
        vec_field_y = zeros(size(M));
    else
        vec_field_x = vx;
        vec_field_y = vy;
    end
    
    M_tilda = imwarp(M, cat(3,vec_field_x,vec_field_y));
    
    [G_fix_x, G_fix_y] = imgradientxy(F,'central');
    [G_fix_mag, ~] = imgradient(G_fix_x, G_fix_y);
    
    iterator = 1;

    mi = mutual_info(F, M);
    
    step_size = 1;
    old_disp_field = cat(3, zeros(size(F)), zeros(size(F)));
    disp_field = cat(3, zeros(size(F)), zeros(size(F)));
    disp_field_diff = Inf;
    ssd = [];

    while iterator<opts.max_iter && disp_field_diff>opts.epsilon
        
        if mod(iterator,100)==0
            ssd = [ssd sum(sum((F - M_tilda).^2))];
            disp(['Iteration number: ' num2str(iterator) ', SSD: ' num2str(ssd(end)) ' and Mutual Info: ' num2str(mi)]);
            step_size = step_size*opts.step;
        end
        
        [G_mov_x, G_mov_y] = imgradientxy(M_tilda, 'central');
        [G_mov_mag, ~] = imgradient(G_mov_x, G_mov_y);
        
        update_field_x = -1 * (((M_tilda - F).*G_fix_x) ./ ...
                ((opts.alpha).^2 * (M_tilda - F).^2 + G_fix_mag.^2) ... 
                + ((M_tilda - F).*G_mov_x) ./ ...
                ((opts.alpha).^2 * (M_tilda - F).^2 + G_mov_mag.^2));
        
        update_field_y = -1 * (((M_tilda - F).*G_fix_y) ./ ...
                ((opts.alpha).^2 * (M_tilda - F).^2 + G_fix_mag.^2) ... 
                + ((M_tilda - F).*G_mov_y) ./ ...
                ((opts.alpha).^2 * (M_tilda - F).^2 + G_mov_mag.^2));
            
        update_field_x(isnan(update_field_x)) = 0;
        update_field_y(isnan(update_field_y)) = 0;
        
        update_field_x = imgaussfilt(update_field_x, opts.sigma_fluid);
        update_field_y = imgaussfilt(update_field_y, opts.sigma_fluid);
        
        [update_field_x,update_field_y] = exponential_mapping(update_field_x,...
                                                            update_field_y);
        if ~opts.compositive
            vec_field_x = vec_field_x + step_size*update_field_x;
            vec_field_y = vec_field_y + step_size*update_field_y;
        else
            [vec_field_x,vec_field_y] = compose_vec_fields(vec_field_x,vec_field_y,...
                step_size*update_field_x,step_size*update_field_y);
        end
        
        vec_field_x = imgaussfilt(vec_field_x, opts.sigma_diff);
        vec_field_y = imgaussfilt(vec_field_y, opts.sigma_diff);
        
        if opts.only_freq; vec_field_x = zeros(size(vec_field_y)); end
        
        disp_field = cat(3, vec_field_x, vec_field_y);
        M_tilda = imwarp(M, disp_field);
        ssd = [ssd, sum(sum((F - M_tilda).^2))];
        iterator = iterator + 1;
        mi = mutual_info(F, M_tilda);
        
        disp_field_diff = sum(sum(sum(abs(old_disp_field - disp_field))));
        old_disp_field = disp_field;
        
        if opts.plot
            subplot(121), imshowpair(F,M_tilda);
            subplot(122), plot(ssd,'r');
            pause(0.001);
        end
    end
    moved_img = M_tilda;
    
    final_MI = mutual_info(F, moved_img);
    final_SSD = sum(sum((F - moved_img).^2));
    
    disp(['Final SSD ' num2str(final_SSD)]);
    disp(['Final mutual info: ' num2str(final_MI)]);
    
%     figure();
%     subplot(221), imshow(F, []), title('Fixed')
%     subplot(222), imshow(M, []), title('Moving')
%     flow = opticalFlow(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)));
%     subplot(223), plot(flow, 'DecimationFactor',[10,5]);
%     subplot(224), imshow(moved_img, [])

%     close all;
%     figure();
%     subplot(131), imshow(F, []), title('Fixed'), subplot(132), imshow(M, []), title('Moving'), ...
%         subplot(133), imshowpair(moved_img, F), title('Moved'), colormap(jet);
end