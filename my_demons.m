function [disp_field,moved_img] = my_demons(fixed_img, moving_img, alpha, sigma_diff, epsilon)

    disp(['Initial SSD ' num2str(sum(sum((fixed_img - moving_img).^2)))]);
    max_iter = 5000;
    vec_field_x = zeros(size(moving_img));
    vec_field_y = zeros(size(moving_img));
    
    [G_fix_x, G_fix_y] = imgradientxy(fixed_img,'central');
    [G_fix_mag, ~] = imgradient(G_fix_x, G_fix_y);
    
    current_moved = moving_img;
    iterator = 1;

    new_mi = mutual_info(fixed_img, moving_img);
    old_mi = new_mi - 100;

    while iterator<max_iter && abs(new_mi - old_mi)>epsilon
        old_mi = new_mi;
        
        if mod(iterator,1000)==0
            disp(['Iteration number: ' num2str(iterator) '   and    ' 'Mutual Info' num2str(new_mi)]);
            disp(['Current mutual info: ' num2str(old_mi - new_mi)]);
        end
        
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
        current_moved = imwarp(moving_img, disp_field);
        iterator = iterator + 1;
        new_mi = mutual_info(fixed_img, current_moved);
    end
    
    disp(['Final mutual info: ' num2str(old_mi - new_mi)]);
    moved_img = current_moved;
    disp(['Initial SSD ' num2str(sum(sum((fixed_img - moved_img).^2)))]);
    
%     subplot(221), imshow(fixed_img, []), title('Fixed')
%     subplot(222), imshow(moving_img, []), title('Moving')
%     flow = opticalFlow(squeeze(disp_field(:,:,1)),squeeze(disp_field(:,:,2)));
%     subplot(223), plot(flow, 'DecimationFactor',[10,5]);
%     subplot(224), imshow(moved_image, [])

    subplot(121), imshowpair(fixed_img, moving_img), title('unregistered')
    subplot(122), imshowpair(fixed_img, moved_img), title('registered')
end