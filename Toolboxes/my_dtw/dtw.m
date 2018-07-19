function [p,q,D] = dtw(M,topology)
    
    [num_rows,num_cols] = size(M);
    D = zeros(num_rows+3,num_cols+3);
    D(:,1:3) = Inf;
    D(1:3,:) = Inf;
    D(1,1) = 0;
    D(2,2) = 0;
    D(3,3) = 0;
    
    phi = zeros(num_rows,num_cols);
    p = num_rows;
    q = num_cols;
    
    if strcmpi(topology,'unconstrained')
        for i = 4:num_rows+3
            for j = 4:num_cols+3
                [local_optimal_cost,tb] = min([D(i-1,j), D(i,j-1), D(i-1,j-1)]);
                D(i,j) = local_optimal_cost + M(i-3,j-3);
                phi(i-3,j-3) = tb;
            end
        end
        
        i = num_rows;
        j = num_cols;
        while i > 1 && j > 1
            tb = phi(i,j);
            if (tb == 1)
                i = i-1;
            elseif (tb == 2)
                j = j-1;
            elseif (tb == 3)
                i = i-1;
                j = j-1;
            else    
                error;
            end
            p = [i,p];
            q = [j,q];
        end
        
    elseif strcmpi(topology,'slopeHalf')
        for i = 4:num_rows+3
            for j = 4:num_cols+3
                [local_optimal_cost,tb] = min([D(i-1,j-2), D(i-2,j-1), D(i-1,j-1)]);
                D(i,j) = local_optimal_cost + M(i-3,j-3);
                phi(i-3,j-3) = tb;
            end
        end
        D(D==Inf) = max(max(D(D~=Inf)));
        
        i = num_rows;
        j = num_cols;
        while i > 1 && j > 1
            tb = phi(i,j);
            if (tb == 1)
                i = i-1;
                j = j-2;
            elseif (tb == 2)
                i = i-2;
                j = j-1;
            elseif (tb == 3)
                i = i-1;
                j = j-1;
            else    
                error;
            end
            p = [i,p];
            q = [j,q];
        end
        
    elseif strcmpi(topology,'slopeThird')
        M = [Inf*ones(num_rows,2) M];
        M = [Inf*ones(2,num_cols+2);M];
        M(1,1) = 0;
        M(2,2) = 0;
        for i = 4:num_rows+3
            for j = 4:num_cols+3
                [local_optimal_cost,tb] = min([D(i-1,j-1)+M(i-1,j-1), D(i-2,j-1)+M(i-2,j-1)+M(i-1,j-1),...
                                          D(i-3,j-1)+M(i-3,j-1)+M(i-3,j-1)+M(i-1,j-1), ...
                                          D(i-1,j-2)+M(i-1,j-2)+M(i-1,j-1) + ...
                                          D(i-1,j-3)+M(i-1,j-3)+M(i-1,j-2)+M(i-1,j-1)]);
                D(i,j) = local_optimal_cost;
                phi(i-3,j-3) = tb;
            end
        end
        D(D==Inf) = max(max(D(D~=Inf)));
        
        i = num_rows;
        j = num_cols;
        while i > 1 && j > 1
            tb = phi(i,j);
            if (tb == 1)
                i = i-1;
                j = j-1;
            elseif (tb == 2)
                i = i-2;
                j = j-1;
            elseif (tb == 3)
                i = i-3;
                j = j-1;
            elseif (tb == 4)
                i = i-1;
                j = j-2;
            elseif (tb == 5)
                i = i-1;
                j = j-3;
            else
                error;
            end
            p = [i,p];
            q = [j,q];
        end
    end
    D = D(4:end, 4:end);
end