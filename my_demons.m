function [disp_field,moved_img,final_SSD,final_MI] = my_demons(fixed_img, moving_img, vx, vy, opts)

%     fixed_img = mat2gray(fixed_img);
%     moving_img = mat2gray(moving_img);
    
    if nargin<3;                        opts                 = struct();     end
    if ~isfield(opts,'alpha');          opts.alpha           = 0.4;          end
    if ~isfield(opts,'sigma_fluid');    opts.sigma_fluid     = 1.0;          end
    if ~isfield(opts,'sigma_diff');     opts.sigma_diff      = 1.0;          end
    if ~isfield(opts,'step');           opts.step            = 1.0;          end
    if ~isfield(opts,'epsilon');        opts.epsilon         = 10;           end
    if ~isfield(opts,'compositive');    opts.compositive     = 0;            end
    if ~isfield(opts,'max_iter');       opts.max_iter        = 1000;         end
    if ~isfield(opts,'plot');           opts.plot            = 0;            end
    

    disp(['Initial SSD ' num2str(sum(sum((fixed_img - moving_img).^2)))]);
    disp(['Initial MI ' num2str(mutual_info(fixed_img, moving_img))]);
    
%     vec_field_x = zeros(size(moving_img));
%     vec_field_y = zeros(size(moving_img));

    vec_field_x = vx;
    vec_field_y = vy;
    current_moved = imwarp(moving_img, cat(3,vec_field_x,vec_field_y));
    
    [G_fix_x, G_fix_y] = imgradientxy(fixed_img,'central');
    [G_fix_mag, ~] = imgradient(G_fix_x, G_fix_y);
    
    iterator = 1;

    mi = mutual_info(fixed_img, moving_img);
    
    step_size = 1;
    old_disp_field = cat(3, zeros(size(fixed_img)), zeros(size(fixed_img)));
    disp_field = cat(3, zeros(size(fixed_img)), zeros(size(fixed_img)));
    disp_field_diff = Inf;
    ssd = [];
    
    while iterator<opts.max_iter && disp_field_diff>opts.epsilon
        
        if mod(iterator,100)==0
            ssd = [ssd sum(sum((fixed_img - current_moved).^2))];
            disp(['Iteration number: ' num2str(iterator) ', SSD: ' num2str(ssd(end)) ' and Mutual Info: ' num2str(mi)]);
            step_size = step_size*opts.step;
        end
        
        [G_mov_x, G_mov_y] = imgradientxy(current_moved, 'central');
        [G_mov_mag, ~] = imgradient(G_mov_x, G_mov_y);
        
        update_field_x = -1 * (((current_moved - fixed_img).*G_fix_x) ./ ...
                ((opts.alpha).^2 * (current_moved - fixed_img).^2 + G_fix_mag.^2) ... 
                + ((current_moved - fixed_img).*G_mov_x) ./ ...
                ((opts.alpha).^2 * (current_moved - fixed_img).^2 + G_mov_mag.^2));
        
        update_field_y = -1 * (((current_moved - fixed_img).*G_fix_y) ./ ...
                ((opts.alpha).^2 * (current_moved - fixed_img).^2 + G_fix_mag.^2) ... 
                + ((current_moved - fixed_img).*G_mov_y) ./ ...
                ((opts.alpha).^2 * (current_moved - fixed_img).^2 + G_mov_mag.^2));
            
        update_field_x(isnan(update_field_x)) = 0;
        update_field_y(isnan(update_field_y)) = 0;
        
        update_field_x = imgaussfilt(update_field_x, opts.sigma_fluid);
        update_field_y = imgaussfilt(update_field_y, opts.sigma_fluid);
        
        if ~opts.compositive
            vec_field_x = vec_field_x + step_size*update_field_x;
            vec_field_y = vec_field_y + step_size*update_field_y;
        else
            [vec_field_x,vec_field_y] = compose_vec_fields(vec_field_x,vec_field_y,...
                step_size*update_field_x,step_size*update_field_y);
        end
        
        vec_field_x = imgaussfilt(vec_field_x, opts.sigma_diff);
        vec_field_y = imgaussfilt(vec_field_y, opts.sigma_diff);
        
        disp_field = cat(3, vec_field_x, vec_field_y);
        current_moved = imwarp(moving_img, disp_field);
        ssd = [ssd, sum(sum((fixed_img - current_moved).^2))];
        iterator = iterator + 1;
        mi = mutual_info(fixed_img, current_moved);
        
        disp_field_diff = sum(sum(sum(abs(old_disp_field - disp_field))));
        old_disp_field = disp_field;
        
        if opts.plot
            subplot(121), imshowpair(fixed_img,current_moved);
            subplot(122), plot(ssd, 'r');
            pause(0.001);
        end
    end
    moved_img = current_moved;
    
    final_MI = mutual_info(fixed_img, moved_img);
    final_SSD = sum(sum((fixed_img - moved_img).^2));
    
    disp(['Final SSD ' num2str(final_SSD)]);
    disp(['Final mutual info: ' num2str(final_MI)]);
    
%     figure();
%     subplot(221), imshow(fixed_img, []), title('Fixed')
%     subplot(222), imshow(moving_img, []), title('Moving')
%     flow = opticalFlow(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)));
%     subplot(223), plot(flow, 'DecimationFactor',[10,5]);
%     subplot(224), imshow(moved_img, [])

%     close all;
%     figure();
%     subplot(131), imshow(fixed_img, []), title('Fixed'), subplot(132), imshow(moving_img, []), title('Moving'), ...
%         subplot(133), imshowpair(moved_img, fixed_img), title('Moved'), colormap(jet);
end