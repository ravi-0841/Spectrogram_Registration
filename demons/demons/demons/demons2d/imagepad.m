%% Pad image
function [I,lim] = imagepad(I,scale)

    if nargin<2; scale = 2; end; % default, pad image twice as big
    
    Ip  = zeros(ceil(size(I)*scale));
    lim = bsxfun(@plus, floor(size(I)*(scale-1)/2), [[1 1];size(I)]); % image limits
    Ip(lim(1):lim(2),lim(3):lim(4)) = I;                              % padded image
    I = Ip;

end