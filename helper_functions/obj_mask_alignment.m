function [new_im_fixed,new_im_movin,moved_mask] = obj_mask_alignment(F,M,psf_size,psf_sigma)
    
    opening_mask = [1;1];
    PSF = fspecial('gaussian', [psf_size psf_size], psf_sigma);
    %% Fixed Image object extraction
    I_trans_F = mat2gray(F);
    I_trans_F = imdiffusefilt(I_trans_F);
    I_decon_F = deconvlucy(I_trans_F, PSF, 10);
    
    level_F = graythresh(I_decon_F);
    I_thresh_F = imbinarize(I_decon_F, level_F);
    
    I_eroded_F = imopen(I_thresh_F, opening_mask);
    
    %% Moving Image object extraction
    I_trans_M = mat2gray(M);
    I_trans_M = imdiffusefilt(I_trans_M);
    I_decon_M = deconvlucy(I_trans_M, PSF, 10);
    
    level_M = graythresh(I_decon_M);
    I_thresh_M = imbinarize(I_decon_M, level_M);
    
    I_eroded_M = imopen(I_thresh_M, opening_mask);

%     optimizer = registration.optimizer.OnePlusOneEvolutionary;
%     metric = registration.metric.MattesMutualInformation;
    [optimizer,metric] = imregconfig('monomodal');
    moved_mask = imregister(mat2gray(I_eroded_M),mat2gray(I_eroded_F),...
                                'translation',optimizer,metric);
%     [~,moved_mask] = imregdemons(mat2gray(I_eroded_M),mat2gray(I_eroded_F),[500 400 100],...
%     'AccumulatedFieldSmoothing',0.6);
%     dsp = my_multires_demons(I_eroded_F, I_eroded_M, 513);
%     moved_mask = imwarp(I_eroded_M, dsp);
    
    level = graythresh(moved_mask);
    moved_mask = imbinarize(moved_mask, level);
    moved_mask = mask_thickening(moved_mask);
    I_eroded_F = mask_thickening(I_eroded_F);
    I_eroded_M = mask_thickening(I_eroded_M);
    
    min_size = 30;
    moved_mask = bwareaopen(moved_mask, min_size, 4);
    I_eroded_F = bwareaopen(I_eroded_F, min_size, 4);
    I_eroded_M = bwareaopen(I_eroded_M, min_size, 4);
    
%     moved_mask(moved_mask>0) = 1;
%     moved_mask = logical(moved_mask);

%%  Scaling Component-wise
    v_fixed = vertical_segmentation(I_eroded_F);
    v_movin = vertical_segmentation(I_eroded_M);
    
    bw_fixed = bwconncomp(v_fixed);
    bw_movin = bwconncomp(v_movin);
    
    assert((bw_movin.NumObjects==bw_fixed.NumObjects),...
            'The utterances may not be same.');
    
    new_im_fixed = zeros(size(I_eroded_F));
    new_im_movin = zeros(size(I_eroded_M));
    
    for obj = 1:bw_movin.NumObjects
        temp_fixed = zeros(bw_fixed.ImageSize);
        temp_movin = zeros(bw_movin.ImageSize);

        temp_fixed(bw_fixed.PixelIdxList{obj}) = 1;
        temp_movin(bw_movin.PixelIdxList{obj}) = 1;
        
        im_fixed_seg = I_eroded_F.*temp_fixed;
        im_movin_seg = I_eroded_M.*temp_movin;
        
        [r_fixed,~] = find(im_fixed_seg==1);
        [r_movin,~] = find(im_movin_seg==1);
        
        ref = max([max(r_fixed), max(r_movin)]);
        
        im_fixed_seg = imresize(F.*I_eroded_F, [513*ref/max(r_fixed), size(im_fixed_seg,2)]);
        im_movin_seg = imresize(M.*I_eroded_M, [513*ref/max(r_movin), size(im_movin_seg,2)]);
        
        im_fixed_seg = im_fixed_seg(1:513,:);
        im_movin_seg = im_movin_seg(1:513,:);
        
        new_im_fixed(temp_fixed==1) = im_fixed_seg(temp_fixed==1);
        new_im_movin(temp_movin==1) = im_movin_seg(temp_movin==1);
    end

    %% Object Mapping
%     moved_cc = bwconncomp(moved_mask,4);
%     fixed_cc = bwconncomp(I_eroded_F,4);
%     movng_cc = bwconncomp(I_eroded_M,4);
%     map_indx = zeros(moved_cc.NumObjects,1);
%     
%     for i = 1:moved_cc.NumObjects
%         if length(moved_cc.PixelIdxList{1,i})<7
%             map_indx(i) = nan;
%             continue;
%         else
%             overlap = zeros(length(fixed_cc.NumObjects),1);
%             for j = 1:fixed_cc.NumObjects
%                 if length(fixed_cc.PixelIdxList{1,j})<7
%                     overlap(j) = 0;
%                     continue;
%                 end
%                 temp_img_fixed = zeros(size(I_eroded_F));
%                 temp_img_moved = zeros(size(moved_mask));
% 
%                 temp_img_moved(moved_cc.PixelIdxList{1,i}) = 1;
%                 temp_img_fixed(fixed_cc.PixelIdxList{1,j}) = 1;
%                 overlap(j) = sum(sum(temp_img_fixed.*temp_img_moved));
%             end
%             if (sum(overlap)==0); map_indx(i) = 0; else; [~, map_indx(i)] ...
%                                                         = max(overlap); end
%         end
%     end
end