function [l, u, W, Bw] = cdpr_ik(P, x)
% inverse kinematics
% x = [x y z phi theta psi]
% l - cable lengths, u - unit vectors, W - structure matrix 6x8

p = x(1:3);
p = p(:);
R = rotm_zyx(x(4), x(5), x(6));

Rb = R*P.B;
Bw = Rb + repmat(p, 1, 8);

Lv = P.A - Bw;              % cable vectors
l = sqrt(sum(Lv.^2, 1))';   % lengths
u = Lv ./ repmat(l', 3, 1);

W = [u; cross(Rb, u, 1)];
end
