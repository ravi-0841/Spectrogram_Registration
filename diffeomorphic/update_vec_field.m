function [ux,uy] = update_vec_field(F,M,init_vx,init_vy,opts)
    M_prime = iminterpolate(M,init_vx,init_vy);
    
    diff = F - M_prime;
    [gy,gx] = imgradientxy(M_prime,'central');
    [grad_mag,~] = imgradient(gx,gy);
    
    scale = diff ./ (grad_mag + diff.^2*opts.sigma_i^2/opts.sigma_x^2);
    scale(grad_mag==0) = 0;
    scale(diff  ==0) = 0;
    ux = gx .* scale;
    uy = gy .* scale;
    
    % Zero non overlapping areas
    ux(F==0)       = 0; uy(F==0)       = 0;
    ux(M_prime==0) = 0; uy(M_prime==0) = 0;
end