function [I_eroded_F,I_eroded_M,moved_mask] = obj_mask_alignment(F,M,psf_size,psf_sigma)
    
    opening_mask = [1;1];
    PSF = fspecial('gaussian', [psf_size psf_size], psf_sigma);
    %% Fixed Image object extraction
    I_trans_F = mat2gray(F);
    I_trans_F = imdiffusefilt(I_trans_F);
    I_decon_F = deconvlucy(I_trans_F, PSF, 10);
    
    level_F = otsuthresh(imhist(I_decon_F));
    I_thresh_F = imbinarize(I_decon_F, level_F);
    
    I_eroded_F = imopen(I_thresh_F, opening_mask);
    
    %% Moving Image object extraction
    I_trans_M = mat2gray(M);
    I_trans_M = imdiffusefilt(I_trans_M);
    I_decon_M = deconvlucy(I_trans_M, PSF, 10);
    
    level_M = otsuthresh(imhist(I_decon_M));
    I_thresh_M = imbinarize(I_decon_M, level_M);
    
    I_eroded_M = imopen(I_thresh_M, opening_mask);
    
    [optimizer,metric] = imregconfig('monomodal');
    moved_mask = imregister(mat2gray(I_eroded_M),mat2gray(I_eroded_F),...
                                'translation',optimizer,metric);
    
    level = graythresh(moved_mask);
    moved_mask = imbinarize(moved_mask, level);
%     moved_mask = mask_thickening(moved_mask);
%     I_eroded_F = mask_thickening(I_eroded_F);
%     I_eroded_M = mask_thickening(I_eroded_M);
%     moved_mask(moved_mask>0) = 1;
%     moved_mask = logical(moved_mask);
                            
    %% Object Mapping
    moved_cc = bwconncomp(moved_mask,4);
    fixed_cc = bwconncomp(I_eroded_F,4);
    movng_cc = bwconncomp(I_eroded_M,4);
    map_indx = zeros(moved_cc.NumObjects,1);
    
    for i = 1:moved_cc.NumObjects
        if length(moved_cc.PixelIdxList{1,i})<7
            map_indx(i) = nan;
            continue;
        else
            overlap = zeros(length(fixed_cc.NumObjects),1);
            for j = 1:fixed_cc.NumObjects
                if length(fixed_cc.PixelIdxList{1,j})<7
                    overlap(j) = 0;
                    continue;
                end
                temp_img_fixed = zeros(size(I_eroded_F));
                temp_img_moved = zeros(size(moved_mask));

                temp_img_moved(moved_cc.PixelIdxList{1,i}) = 1;
                temp_img_fixed(fixed_cc.PixelIdxList{1,j}) = 1;
                overlap(j) = sum(sum(temp_img_fixed.*temp_img_moved));
            end
            if (sum(overlap)==0); map_indx(i) = 0; else; [~, map_indx(i)] ...
                                                        = max(overlap); end
        end
    end
end