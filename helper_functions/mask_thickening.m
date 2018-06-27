function out = mask_thickening(im)
    out = im;
    input_conn = bwconncomp(im,4);
    ref_num = input_conn.NumObjects;
    red_img = im(2:200,2:end-2);
    [r,c] = find(red_img==0);
    r = r + 1;
    c = c + 1;
    for idx = 1:length(r)
        temp = im;
%         if (im(r(idx)-1,c(idx))==1) || (im(r(idx),c(idx)-1)==1) || ...
%                 (im(r(idx),c(idx)+1)==1) || (im(r(idx)+1,c(idx))==1)
        if (im(r(idx)-1,c(idx))==1) || (im(r(idx),c(idx)-1)==1)
            temp(r(idx),c(idx)) = 1;
            temp_conn = bwconncomp(temp,4);
            if temp_conn.NumObjects < ref_num
                continue;
            else
                out(r(idx), c(idx)) = 1;
            end
        end
    end
end