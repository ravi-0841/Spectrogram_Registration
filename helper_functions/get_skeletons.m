function shapes = get_skeletons(bw_img)
    shapes = [];
    cc = bwconncomp(bw_img);
    for i = 1:cc.NumObjects
        z = zeros(cc.ImageSize);
        z(cc.PixelIdxList{i}) = 1;
        skeleton = bwskel(logical(z));
        [r,c] = find(skeleton==1);
        knots = [c'-mean(c);r'-mean(r)];
        num_pts = length(r);
        original_scale = 1:num_pts;
        finer_scale = linspace(1,num_pts,50);
        splineXY = spline(original_scale, knots, finer_scale);
        shapes = [shapes [splineXY(1,:)' ; splineXY(2,:)']];
    end
end