function [r, w] = LDL_solve(X, y)
% Input:
%   - X \in R^{m x n}
%   - y \in R^n 
% Output:
%   - r \in R^m such that r + X * w = y  and X^T * r = 0
%   - w \in R^n is the solution of the least square problem \min_{w} || Xw - y ||
[m, n] = size(X);
% Compute A
A = zeros(m + n, m + n);
A(1:m, 1:m) = eye(m, m);
A(m+1:end, 1:m) = X';
A(1:m, m+1:end) = X;
% Compute b
b = zeros(m + n, 1);
b(1:m) = y;
% Compute the LDL factorization of A
[L, D, perm, pivot] = LDL_factorization(A);
% Since L is lower triangular, we use a solver designed for lower triangular matrices
opts.LT = true;
z = linsolve(L, b(perm), opts);
q = solve_block_diagonal(D, z, pivot);
% Since L^T is upper triangular, we use a solver designed for upper triangular matrices
opts.LT = false;
opts.UT = true;
t = linsolve(L', q, opts);
perm_t = compute_P_transpose(perm);   
x = t(perm_t)';
r = x(1:m);
w = x(m+1:end);
end

