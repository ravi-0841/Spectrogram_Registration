function out = vertical_segmentation(bin_image)
    out = zeros(size(bin_image));
    [~,c] = find(bin_image==1);
    non_zero_cols = sum(bin_image(:,c));
    for i = 2:length(non_zero_cols)-1
        if non_zero_cols(i)<4 && non_zero_cols(i-1)<4 ...
                && non_zero_cols(i+1)<4
            non_zero_cols(i) = 0;
        end
    end
%     idx = find(non_zero_cols>1);
    out(:,c(non_zero_cols>0)) = 1;
end