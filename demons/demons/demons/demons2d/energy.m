%% Get energy
function e = energy(F,M,sx,sy,sigma_i,sigma_x)

    % Intensity difference
    Mp     = iminterpolate(M,sx,sy);
    diff2  = (F-Mp).^2;
    area   = size(M,1)*size(M,2);
    
    % Transformation Gradient
    jac = jacobian(sx,sy);
    
    % Three energy components
    e_sim  = sum(diff2(:)) / area;
    %e_dist = sum((cx(:)-sx(:)).^2 + (cy(:)-sy(:)).^2) / area;
    e_reg = sum(jac(:).^2) / area;
    
    % Total energy
    e      = e_sim + (sigma_i^2/sigma_x^2) * e_reg;

end