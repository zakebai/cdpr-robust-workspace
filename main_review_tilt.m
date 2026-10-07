%% R5: are the 19 tilt scenarios representative?
% reference robot, z = 3 m, payload 350-650 kg + wind 25 m/s
% compare: 19 scenarios (paper) vs 19 scenarios + 100 random tilt sets (each mast, each axis uniform in +-1 deg)
clear; warning('off','all');
P = cdpr_params(1,'crossed');
Fw = 0.5*1.225*1.2*0.96*25^2;
WC = cdpr_wind_vertices(P, [350 650], Fw);
g = deg2rad(1); Pt = {P};
for j = 1:4, for d = 1:2, for sg = [1 -1]
    G = zeros(4,2); G(j,d) = sg*g; Q = P; Q.A = cdpr_tilt_all(P, G); Pt{end+1} = Q;
end, end, end
Q = P; Q.A = cdpr_tilt_all(P,  g*ones(4,2)); Pt{end+1} = Q;
Q = P; Q.A = cdpr_tilt_all(P, -g*ones(4,2)); Pt{end+1} = Q;
rand('seed', 7); Pr = {};
for k = 1:100
    Q = P; Q.A = cdpr_tilt_all(P, g*(2*rand(4,2)-1)); Pr{end+1} = Q;
end
nx = 30; dx = 2*P.L/nx; xs = -P.L+dx/2:dx:P.L; ny = 20; dy = 2*P.W/ny; ys = -P.W+dy/2:dy:P.W;
n19 = 0; nboth = 0; t0 = tic;
for i = 1:ny
  for j = 1:nx
    x = [xs(j); ys(i); 3];
    if cdpr_robust_ok(Pt, x, WC)
        n19 = n19 + 1;
        nboth = nboth + cdpr_robust_ok(Pr, x, WC);
    end
  end
end
fprintf('R5 z=3, %d cells: 19 scenarios %.1f %% | 19 + 100 random %.1f %% | %.0f s\n', nx*ny, 100*n19/(nx*ny), 100*nboth/(nx*ny), toc(t0));
