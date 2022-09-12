function [m, row, col] = max_abs_diagonal(A)
% Input: A \in R^{m x n}
% Output:
%   - m = maximum absolute value of diag(A)
%   - row = row of the element having maximum absolute value of diag(A)
%   - col = column of the element having maximum absolute value of diag(A)
[m, row] = max(abs(diag(A)));
col = row;
end

