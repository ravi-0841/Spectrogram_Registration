function [disp_field,moved_image] = demons_registration(fixed_img, moving_img, alpha, sigma_diff, max_iter)
    vec_field_x = zeros(size(moving_img));
    vec_field_y = zeros(size(moving_img));
    
    [G_fix_x, G_fix_y] = imgradientxy(fixed_img,'central');
    [G_fix_mag, ~] = imgradient(G_fix_x, G_fix_y);
    
    current_moved = moving_img;
    
    for i = 1:max_iter
%         disp(['Iteration number: ' num2str(i)]);
        [G_mov_x, G_mov_y] = imgradientxy(current_moved, 'central');
        [G_mov_mag, ~] = imgradient(G_mov_x, G_mov_y);
        
        vec_field_x = vec_field_x + ((current_moved - fixed_img).*G_fix_x) ./ ...
                (alpha.^2 * (current_moved - fixed_img).^2 + G_fix_mag.^2) ... 
                + ((current_moved - fixed_img).*G_mov_x) ./ ...
                (alpha.^2 * (current_moved - fixed_img).^2 + G_mov_mag.^2);
        
        vec_field_y = vec_field_y + ((current_moved - fixed_img).*G_fix_y) ./ ...
                (alpha.^2 * (current_moved - fixed_img).^2 + G_fix_mag.^2) ... 
                + ((current_moved - fixed_img).*G_mov_y) ./ ...
                (alpha.^2 * (current_moved - fixed_img).^2 + G_mov_mag.^2);
        
        vec_field_x = imgaussfilt(vec_field_x, sigma_diff);
        vec_field_y = imgaussfilt(vec_field_y, sigma_diff);
        
        disp_field = cat(3, vec_field_x, vec_field_y);
        current_moved = imwarp(current_moved, disp_field);
        
    end
    moved_image = current_moved;
    
    flow = opticalFlow(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)));
    subplot(121), plot(flow, 'DecimationFactor',[10,5]);
    subplot(122), imshow(moved_image, [])
end