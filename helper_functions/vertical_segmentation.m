function out = vertical_segmentation(bin_image)
    out = zeros(size(bin_image));
    [~,c] = find(bin_image==1);
    out(:,c) = 1;
end