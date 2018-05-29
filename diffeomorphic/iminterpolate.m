%% Interpolate image
function I = iminterpolate(I,vx,vy)

    [x,y] = ndgrid(0:(size(I,1)-1), 0:(size(I,2)-1));
    x_prime = x + vx;
    y_prime = y + vy;

    I = interpn(x,y,I,x_prime,y_prime,'linear',0);
end