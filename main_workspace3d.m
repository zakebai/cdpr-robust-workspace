%% 3D wrench feasible workspace with orientation sets (reference robot)
% grid: cell centres over the mast footprint 12 x 8 m, z = 0..5.5 m
% results: volume and % for 4 orientation sets, crossed and parallel

clear; clc; close all;
warning('off', 'all');

nx = 40;
topo = {'crossed', 'parallel'};
r1 = 0.1; r3 = 0.3;           % roll/pitch and yaw ranges, rad
S{1} = [0;0;0];
S{2} = [0 0 0; 0 0 0; 0 r3 -r3];
S{3} = [0 r1 -r1 0 0; 0 0 0 r1 -r1; 0 0 0 0 0];
S{4} = [S{2}, S{3}(:,2:end)];
sname = {'zero', 'yaw +-0.3', 'roll/pitch +-0.1', 'all'};

P = cdpr_params(1, 'crossed');
dx = 2*P.L/nx;
xs = -P.L + dx/2 : dx : P.L;
ny = round(nx*P.W/P.L);
dy = 2*P.W/ny;
ys = -P.W + dy/2 : dy : P.W;
dz = 0.5;
zs = dz/2 : dz : 5.5;
Vbox = (2*P.L)*(2*P.W)*5.5;

V = zeros(2, 4);
WS = cell(2, 4);
for t = 1:2
    P = cdpr_params(1, topo{t});
    w = cdpr_wrench(P);
    for s = 1:4
        F = false(length(ys), length(xs), length(zs));
        for k = 1:length(zs)
            for i = 1:length(ys)
                for j = 1:length(xs)
                    F(i,j,k) = cdpr_robust_ok({P}, [xs(j); ys(i); zs(k)], w, S{s});
                end
            end
            fprintf('.');
        end
        fprintf('\n');
        WS{t,s} = F;
        V(t,s) = sum(F(:))*dx*dy*dz;
        fprintf('%-8s %-18s volume %6.1f m3 (%5.1f %% of 12x8x5.5 box)\n', ...
            topo{t}, sname{s}, V(t,s), 100*V(t,s)/Vbox);
    end
end

% area of each slice (crossed, zero orientation and all orientations)
Az = zeros(length(zs), 2);
for k = 1:length(zs)
    Az(k,1) = 100*mean(mean(WS{1,1}(:,:,k)));
    Az(k,2) = 100*mean(mean(WS{1,4}(:,:,k)));
end
disp('z, m | % of footprint (zero orientation) | % (all orientations)');
disp([zs' Az]);

% convergence of the footprint % at z = 3 m (crossed, zero orientation)
P = cdpr_params(1, 'crossed');
w = cdpr_wrench(P);
conv = [];
for n = [20 40 80 160]
    d = 2*P.L/n;
    m = round(n*P.W/P.L);
    xg = -P.L + d/2 : d : P.L;
    yg = -P.W + P.W/m : 2*P.W/m : P.W;
    c = 0;
    for i = 1:length(yg)
        for j = 1:length(xg)
            [~, ~, W] = cdpr_ik(P, [xg(j); yg(i); 3; 0; 0; 0]);
            c = c + cdpr_feasible(W, w, P.tmin, P.tmax);
        end
    end
    conv = [conv; n, m, 100*c/(n*m)];
    fprintf('cells %d x %d: %.1f %% of footprint at z = 3 m\n', n, m, 100*c/(n*m));
end

save results_ws3d.mat V WS xs ys zs Az sname Vbox conv

%% Figure: slices of the crossed robot
figure('Position', [100 100 900 350]);
for k = 1:3
    subplot(1,3,k);
    kk = [2 6 10];
    imagesc(xs, ys, WS{1,1}(:,:,kk(k)) + WS{1,4}(:,:,kk(k))); axis equal tight; set(gca, 'YDir', 'normal');
    title(sprintf('z = %.2f m', zs(kk(k))));
end
print('-dpng', '-r200', 'fig_ws3d.png');
