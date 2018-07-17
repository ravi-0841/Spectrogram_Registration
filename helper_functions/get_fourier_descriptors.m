function shapes = get_fourier_descriptors(img, num_coeffs)
    bwcc    = bwconncomp(img);
    shapes  = zeros(num_coeffs*4, bwcc.NumObjects);
    for obj = 1:bwcc.NumObjects
        z                           = zeros(bwcc.ImageSize);
        z(bwcc.PixelIdxList{obj})   = 1;
        [a,b,c,d,~]                 = ellipticFourierDescriptor(z,num_coeffs,1,0);
        shapes(:,obj)               = [a;b;c;d];
    end
end