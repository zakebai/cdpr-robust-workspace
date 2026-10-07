function P = cdpr_params(s, topo)
% parameters of 8 cable robot
% s - scale (1 = big model, 0.25 = lab model)
% topo - 'crossed' or 'parallel'

if nargin < 1, s = 1; end
if nargin < 2, topo = 'crossed'; end

P.scale = s;
P.topology = topo;

% frame size (half length, half width, height)
P.L = 6*s;
P.W = 4*s;
P.H = 6*s;
P.d = 0.5*s;   % distance between 2 pulleys on one mast

% platform size (half)
P.pl = 0.8*s;
P.pw = 0.5*s;
P.ph = 0.3*s;

% mass and tension limits
P.m = 500*s^3;
P.g = [0; 0; -9.81];
P.tmin = 0.02*P.m*9.81;
if s == 1
    P.tmin = 100;
end
P.tmax = 15000*s^3;
P.r_pulley = 0.1*s;

% 4 masts
towers = [1 1; -1 1; -1 -1; 1 -1];

A = zeros(3,8);
B = zeros(3,8);
for j = 1:4
    sx = towers(j,1);
    sy = towers(j,2);
    if strcmp(topo, 'crossed')
        A(:,2*j-1) = [sx*(P.L-P.d); sy*P.W; P.H];
        A(:,2*j)   = [sx*P.L; sy*(P.W-P.d); P.H];
        B(:,2*j-1) = [sx*P.pl; -sy*P.pw; P.ph];
        B(:,2*j)   = [-sx*P.pl; sy*P.pw; -P.ph];
    else
        A(:,2*j-1) = [sx*P.L; sy*P.W+P.d; P.H];
        A(:,2*j)   = [sx*P.L; sy*P.W-P.d; P.H];
        B(:,2*j-1) = [sx*P.pl; sy*P.pw; P.ph];
        B(:,2*j)   = [sx*P.pl; sy*P.pw; -P.ph];
    end
end

P.A = A;   % pulley points
P.B = B;   % points on platform
P.towers = towers;
end
