function [vx, vy] = exponential_mapping(vx, vy)

    vec_mag = vx.^2 + vy.^2;
    m = sqrt(max(vec_mag(:)));
    n = ceil(log2(m/0.5));
    n = max(n,0);

    vx = vx * 2^-n;
    vy = vy * 2^-n;

    for i=1:n
        [vx,vy] = compose(vx,vy, vx,vy);
    end

end