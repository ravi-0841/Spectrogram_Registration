function M = pair_sim(A,B,metric)

EA = sqrt(sum(A.^2));
EB = sqrt(sum(B.^2));

if strcmpi(metric,'Cosine')
    M = (A'*B)./(EA'*EB);
else
    ncA = size(A,2);
    ncB = size(B,2);
    M = zeros(ncA, ncB);
    for i = 1:ncA
     for j = 1:ncB
       % Euclidean Distance
       M(i,j) = (sum((A(:,i) - B(:,j)).^2));
     end
    end
end