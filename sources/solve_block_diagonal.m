function [w] = solve_block_diagonal(D, z, pivot)
% Input:
%   - D \in R^{n x n} block-diagonal
%   - z \in R^n

% Output:
%   - w \in R^n is the solution of the linear system D * w = z solved 
%   by exploiting the fact that D is block-diagonal.

n = size(D, 1);
w = zeros(n);

k = 1;
while(k <= n)
    if(pivot(k) == 1)
        w(k) = z(k) / D(k,k);
        k = k + 1;
    else
        w(k:k+1) = linsolve(D(k:k+1, k:k+1), z(k:k+1));
        k = k + 2;
    end

end
end

