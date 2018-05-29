function [disp_field,moved_img,final_SSD,final_MI] = my_demons(fixed_img, moving_img, opts)

    if nargin<3;                        opts                 = struct();     end
    if ~isfield(opts,'alpha');          opts.alpha           = 0.4;          end
    if ~isfield(opts,'sigma_fluid');    opts.sigma_fluid     = 1.0;          end
    if ~isfield(opts,'sigma_diff');     opts.sigma_diff      = 1.0;          end
    if ~isfield(opts,'step');           opts.step            = 1.0;          end
    if ~isfield(opts,'epsilon');        opts.epsilon         = 0.0001;       end
    

    disp(['Initial SSD ' num2str(sum(sum((fixed_img - moving_img).^2)))]);
    disp(['Initial MI ' num2str(mutual_info(fixed_img, moving_img))]);
    
    max_iter = 2000;
    vec_field_x = zeros(size(moving_img));
    vec_field_y = zeros(size(moving_img));
    
    [G_fix_x, G_fix_y] = imgradientxy(fixed_img,'central');
    [G_fix_mag, ~] = imgradient(G_fix_x, G_fix_y);
    
    current_moved = moving_img;
    iterator = 1;

    new_mi = mutual_info(fixed_img, moving_img);
    old_mi = new_mi - 100;
    
    step_size = 1;
    old_disp_field = cat(3, zeros(size(fixed_img)), zeros(size(fixed_img)));
    disp_field = cat(3, zeros(size(fixed_img)), zeros(size(fixed_img)));
    disp_field_diff = Inf;

    while iterator<max_iter && disp_field_diff>opts.epsilon
        old_mi = new_mi;
        
        if mod(iterator,500)==0
            ssd = sum(sum((fixed_img - current_moved).^2));
            disp(['Iteration number: ' num2str(iterator) ', SSD: ' num2str(ssd) ' and Mutual Info: ' num2str(new_mi)]);
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
        
        vec_field_x = vec_field_x + step_size*update_field_x;
        vec_field_y = vec_field_y + step_size*update_field_y;
        
        vec_field_x = imgaussfilt(vec_field_x, opts.sigma_diff);
        vec_field_y = imgaussfilt(vec_field_y, opts.sigma_diff);
        
        disp_field = cat(3, vec_field_x, vec_field_y);
        current_moved = imwarp(moving_img, disp_field);
%         current_moved = iminterpolate(moving_img,vec_field_x,vec_field_y);
        iterator = iterator + 1;
        new_mi = mutual_info(fixed_img, current_moved);
        
        disp_field_diff = sum(sum(sum(abs(old_disp_field - disp_field))));
        old_disp_field = disp_field;
        imshow(current_moved, []), colormap(jet)
        pause(0.01);
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

%     figure()
%     subplot(121), imshowpair(moving_img, fixed_img), title('unregistered')
%     subplot(122), imshowpair(moved_img, fixed_img), title('registered')

    figure()
    subplot(131), imshow(fixed_img, []), title('Fixed'), subplot(132), imshow(moving_img, []), title('Moving'), ...
        subplot(133), imshow(moved_img, []), title('Moved'), colormap(jet);
end