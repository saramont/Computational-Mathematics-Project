function [Q1,R1] = QR_factorization(X)
% Computes the thin-qr factorization of a tall-thin matrix X.
%Input: X \in R^{m x n} (m >> n)
%Output: 
%   - R1 \in R^{n x n}
%   - Q1 \in R^{m x n}

[m, n] = size(X);
U = zeros(m, n);
for j = 1:n
    [U(j:end, j), X(j,j)] = householder_vector(X(j:end, j));
    X(j+1:end,j) = 0;
    X(j:end,j+1:end) = X(j:end,j+1:end) - 2 * U(j:end, j) * (U(j:end, j)' * X(j:end,j+1:end));
end
R1 = X(1:n, 1:end);
% Compute Q1 using the Householder reflectors
Q1 = eye(m,n);
for j = n:-1:1
    Q1 = Q1 - 2 * U(1:end, j) * (U(1:end, j)' * Q1);
end
end


