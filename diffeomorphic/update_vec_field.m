function [ux,uy] = update_vec_field(F,M,init_vx,init_vy,opts)
    M_prime = iminterpolate(M,init_vx,init_vy);
%     M_prime = imwarp(M, cat(3, init_vx, init_vy));
    
    diff = F - M_prime;
    [gy,gx] = imgradientxy(M_prime,'central');
    [grad_mag,~] = imgradient(gx,gy);
    
    scale = (diff ./ (grad_mag + diff.^2*opts.sigma_i^2/opts.sigma_x^2));
    scale(isnan(scale)) = 0;
    
    ux = gx .* scale;
    uy = gy .* scale;
end