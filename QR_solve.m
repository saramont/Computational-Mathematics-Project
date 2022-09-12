function [w] = QR_solve(X, y)
% Inputs:
%   - X \in R^{m x n}
%   - y \in R^m is a random vector
% Output: w \in R^n is the solution of the least square problem \min_{w} || Xw - y ||
[Q1, R1] = QR_factorization(X);
c = Q1' * y;
% Since R1 is upper triangular, we use a solver designed for upper triangular matrices
opts.UT = true;
w = linsolve(R1,c, opts);
end

