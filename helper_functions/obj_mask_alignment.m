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
end