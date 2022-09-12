function [u, s] = householder_vector(x)
% Input: vector x. 
% Output: 
%   - u = housolder vector of x
%   - s = 2-norm of x
s = norm(x);
if x(1) >= 0, s = -s; end
v = x;
v(1) = v(1) - s;
u = v / norm(v);

