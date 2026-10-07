%% Combined worst case: mass + wind + mast tilt together (reference robot)
% A nominal: 500 kg, no wind, masts straight
% B mass 350..650 kg
% C mass 350..650 kg + wind 25 m/s from any horizontal direction
% D case C + mast tilt up to 1 deg (19 tilt cases)
% result: % of the 12 x 8 m footprint at z = 1..5 m

clear; clc; close all;
warning('off', 'all');

P = cdpr_params(1, 'crossed');
nx = 60;
dx = 2*P.L/nx;
xs = -P.L + dx/2 : dx : P.L;
ny = round(nx*P.W/P.L);
dy = 2*P.W/ny;
ys = -P.W + dy/2 : dy : P.W;
zs = 1:5;

rho = 1.225; Cd = 1.2; Ar = 0.96; v = 25;
Fw = 0.5*rho*Cd*Ar*v^2;
fprintf('wind force %.0f N\n', Fw);

WA = cdpr_wind_vertices(P, 500, 0);
WB = cdpr_wind_vertices(P, [350 650], 0);
WC = cdpr_wind_vertices(P, [350 650], Fw);

% tilt cases: each mast 1 deg in +-X or +-Y, and all masts in/out
g = deg2rad(1);
Pt = {P};
for j = 1:4
    for d = 1:2
        for sg = [1 -1]
            G = zeros(4,2); G(j,d) = sg*g;
            Q = P; Q.A = cdpr_tilt_all(P, G); Pt{end+1} = Q;
        end
    end
end
Q = P; Q.A = cdpr_tilt_all(P,  g*ones(4,2)); Pt{end+1} = Q;
Q = P; Q.A = cdpr_tilt_all(P, -g*ones(4,2)); Pt{end+1} = Q;
fprintf('%d geometries in case D\n', numel(Pt));

R = zeros(length(zs), 4);
MAP = cell(length(zs), 1);
for k = 1:length(zs)
    M = zeros(length(ys), length(xs));     % 0 = out, 1..4 = last case passed
    for i = 1:length(ys)
        for j = 1:length(xs)
            x = [xs(j); ys(i); zs(k)];
            if ~cdpr_robust_ok({P}, x, WA), continue; end
            M(i,j) = 1;
            if ~cdpr_robust_ok({P}, x, WB), continue; end
            M(i,j) = 2;
            if ~cdpr_robust_ok({P}, x, WC), continue; end
            M(i,j) = 3;
            if ~cdpr_robust_ok(Pt, x, WC), continue; end
            M(i,j) = 4;
        end
    end
    MAP{k} = M;
    for c = 1:4
        R(k,c) = 100*mean(M(:) >= c);
    end
    fprintf('z = %d m: A %.1f | B %.1f | C %.1f | D %.1f %%\n', zs(k), R(k,:));
end

save results_worstcase.mat R MAP xs ys zs Fw

figure('Position', [100 100 600 400]);
imagesc(xs, ys, MAP{3}); axis equal tight; set(gca, 'YDir', 'normal'); colorbar;
title('z = 3 m: 1 nominal, 2 +mass, 3 +wind, 4 +tilt');
print('-dpng', '-r200', 'fig_worstcase.png');
