function m = mutual_info(A, B)
    A = uint8(255*(mat2gray(A)));
    B = uint8(255*(mat2gray(B)));
    
    hA = imhist(A);
    hB = imhist(B);
    
    hA = hA./sum(hA);
    hB = hB./sum(hB);
    
    hA(hA==0) = 1;
    hB(hB==0) = 1;
    
    entA = -1 * sum(hA.*log2(hA));
    entB = -1 * sum(hB.*log2(hB));
    
    M = zeros(256,256);
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            try
                M(A(i,j)+1,B(i,j)+1) = M(A(i,j)+1,B(i,j)+1) + 1;
            catch error
                disp([A(i,j)+1 B(i,j)+1]);
            end
        end
    end
    M = M./sum(sum(M));
    M(M==0) = 1;
    jointAB = -1 * sum(sum(log2(M).*M));
    m = entA + entB - jointAB;
end