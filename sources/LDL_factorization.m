function [L, A, perm, pivot] = LDL_factorization(A)
% Input: A \in R^{n x n} symmetric.
% Output:
%   - L \in R^{n x n} lower triangular matrix (with ones in the diagonal).
%   - A \in R^{n x n} block-diagonal matrix.
%   - perm = vector reppresenting the permutation matrix P.
%   - pivot = vector used to store the pivot block's size at each step
% It computes the LDL factorization of P*A*P' = L*D*L' (the algorithm modifies A so that in the end A = D).
n = size(A,1);
L = eye(n);
pivot = zeros(1, n);
perm = 1:n;
alpha = (1 + sqrt(17)) / 8;
k = 1;
while(k <= (n-1))
    % Compute mu_0 and its indexes
    [mu_0, r_0, c_0] = max_abs_matrix(A(k:end, k:end));
    r_0 = r_0 + k - 1;
    c_0 = c_0 + k - 1;
    % Compute mu_1 and its indexes
    [mu_1, r_1, c_1] = max_abs_diagonal(A(k:end, k:end)); % Note: r_1 = c_1 
    r_1 = r_1 + k - 1;
    c_1 = c_1 + k - 1;
    if(mu_0 <= 1e-16 && mu_1 <= 1e-16)
        error('The input matrix is singular: cannot compute the LDL factorization.')
    end
    if(mu_1 >= (alpha * mu_0))
        % 1x1 pivoting 
        % Exchange row k with row r_1 in A
        A([k, r_1], 1:end) = A([r_1, k], 1:end);
        % Exchange column k with column c_1 in A
        A(1:end, [k, c_1]) = A(1:end, [c_1, k]);
        
        % Update the permutation vector
        perm([k, r_1]) = perm([r_1, k]);
        
        % Exchange row k with row r_1 in L
        L([k, r_1], 1:k-1) = L([r_1, k], 1:k-1);
        % Update L
        L(k+1:end, k) = A(k+1:end, k) / A(k, k);
        
        % Update the submatrix
        A(k+1:end, k+1:end) = A(k+1:end, k+1:end) - L(k+1:end, k) * A(k, k+1:end);
        A(k+1:end, k) = 0;
        A(k, k+1:end) = 0;
        
        pivot(k) = 1;
        k = k + 1;
    else
        % 2x2 pivoting
        % Exchange row k with row c_0 in A
        A([k, c_0], 1:end) = A([c_0, k], 1:end);
        % Exchange column k with column c_0 in A
        A(1:end, [k, c_0]) = A(1:end, [c_0, k]);
        
        % Update the permutation vector
        perm([k, c_0]) = perm([c_0, k]);
        
        % Exchange row k+1 with row r_0 in A
        A([k+1, r_0], 1:end) = A([r_0, k+1], 1:end);
        % Exchange column k+1 with column r_0 in A
        A(1:end, [k+1, r_0]) = A(1:end, [r_0, k+1]);     
        
        % Update the permutation vector
        perm([k+1, r_0]) = perm([r_0, k+1]);     
        
        % Exchange row k with row c_0 in L
        L([k, c_0], 1:k-1) = L([c_0, k], 1:k-1);
        % Exchange row k+1 with row r_0 in L
        L([k+1, r_0], 1:k-1) = L([r_0, k+1], 1:k-1);
        % Update L
        L(k+2:end, k:k+1) = A(k+2:end, k:k+1) / A(k:k+1, k:k+1);
        
        % Update the submatrix
        A(k+2:end, k+2:end) = A(k+2:end, k+2:end) - L(k+2:end, k:k+1) * A(k:k+1, k+2:end);
        
        A(k+2:end, k:k+1) = 0;
        A(k:k+1, k+2:end) = 0;
        
        pivot(k) = 2;
        k = k + 2;
    end
end
% Update the pivot information in position n
% pivot(n-1) = 2 -> pivot(n) = 0
% pivot(n-1) = 1 -> pivot(n) = 1
% pivot(n-1) = 0 -> pivot(n) = 1
if(pivot(n-1) ~= 2)
    pivot(n) = 1;
end
end

