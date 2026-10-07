function w = cdpr_wrench(P, a, F, M)
% external wrench, W*tau = w
% a - acceleration, F - force, M - moment

if nargin < 2, a = [0;0;0]; end
if nargin < 3, F = [0;0;0]; end
if nargin < 4, M = [0;0;0]; end

w = [P.m*(a(:) - P.g) - F(:); -M(:)];
end
