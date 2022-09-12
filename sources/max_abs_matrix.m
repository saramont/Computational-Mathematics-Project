function [value, row, col] = max_abs_matrix(A)
% Input: A \in R^{m x n}
% Output:
%   - m = maximum absolute value of A
%   - row = row of the element having maximum absolute value of A
%   - col = column of the element having maximum absolute value of A
[values, rows] = max(abs(A));
[value, col] = max(values);
row = rows(col);
end

