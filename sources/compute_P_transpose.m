function [perm_t] = compute_P_transpose(perm)
% Input: perm = a vector rappresenting a permutation matrix P
% Output: perm_t = the vector rappresenting the permutation matrix P'
[~, n] = size(perm);
perm_t = zeros(1, n);
for i = 1:n
    for j = 1:n
        if(i == perm(j))
           perm_t(i) = j;
        end
    end
end      
end

