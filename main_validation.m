%% Validation of the method
% 1. kinematics check on the real CoGiRo geometry
% 2. two different algorithms for workspace give the same answer
% 3. CoGiRo static capacity vs documented payload
% 4. grid convergence

clear; clc; close all;
warning('off', 'all');
rng(1);

%% 1. CoGiRo kinematics: inverse -> forward
Pc = cogiro_params();
err_fk = 0; err_J = 0;
for k = 1:100
    x = [-5 + 10*rand; -3.5 + 7*rand; 0.5 + 3*rand; 0.3*(rand(3,1) - 0.5)];
    l = cdpr_ik(Pc, x);
    xe = cdpr_fk(Pc, l, x + [0.05*randn(3,1); 0.01*randn(3,1)]);
    err_fk = max(err_fk, norm(xe - x));

    % jacobian check: dl = -W'*[dp; dw] for small motion
    % (dp - small translation, dw - small rotation vector)
    [~, ~, W] = cdpr_ik(Pc, x);
    h = 1e-6;
    dp = h*[1; 0.5; -0.3];
    dw = h*[0.3; -0.2; 0.5];
    R = rotm_zyx(x(4), x(5), x(6));
    S = [0 -dw(3) dw(2); dw(3) 0 -dw(1); -dw(2) dw(1) 0];
    R2 = expm(S)*R;                               % rotated platform
    Lv = Pc.A - repmat(x(1:3) + dp, 1, 8) - R2*Pc.B;
    l2 = sqrt(sum(Lv.^2, 1))';
    err_J = max(err_J, norm((l2 - l) - (-W'*[dp; dw]))/h);
end
fprintf('CoGiRo: FK error max %.1e, jacobian error %.1e\n', err_fk, err_J);

%% 2. LP vs capacity margin (hyperplane shifting)
P = cdpr_params(1, 'crossed');
Npt = 2000;
agree = 0; nfeas = 0; minabs = inf;
for k = 1:Npt
    if k <= Npt/2
        Q = P;  name = 'ref';
        x = [-6 + 12*rand; -4 + 8*rand; 0.5 + 5*rand; 0.4*(rand(3,1) - 0.5)];
        w = cdpr_wrench(Q, [0;0;0], 441*[cos(2*pi*rand); sin(2*pi*rand); 0]);
    else
        Q = Pc;
        x = [-7 + 14*rand; -5 + 10*rand; 0.5 + 4*rand; 0.4*(rand(3,1) - 0.5)];
        R = rotm_zyx(x(4), x(5), x(6));
        mt = Q.m + 600*rand;
        w = mt*[-Q.g; -cross(R*Q.com, Q.g)];
    end
    [~, ~, W] = cdpr_ik(Q, x);
    a = cdpr_feasible(W, w, Q.tmin, Q.tmax, 1);    % LP only
    s = cdpr_capacity(W, w, Q.tmin, Q.tmax);
    b = (s >= 0);
    agree = agree + (a == b);
    nfeas = nfeas + b;
    minabs = min(minabs, abs(s));
end
fprintf('LP and capacity margin agree in %d of %d poses (%d feasible), smallest |margin| %.2f N\n', agree, Npt, nfeas, minabs);

%% 3. CoGiRo: static payload capacity (no safety factor, no dynamics)
zz = [1 2 3];
xs = linspace(-6, 6, 25);
ys = linspace(-4, 4, 17);
pay = zeros(length(ys), length(xs), length(zz));
for kz = 1:length(zz)
    for i = 1:length(ys)
        for j = 1:length(xs)
            pay(i,j,kz) = cdpr_maxload(Pc, [xs(j); ys(i); zz(kz); 0; 0; 0]) - Pc.m;
        end
    end
end
for kz = 1:length(zz)
    pz = pay(:,:,kz);
    pc = cdpr_maxload(Pc, [0; 0; zz(kz); 0; 0; 0]) - Pc.m;
    fprintf('CoGiRo z = %d m: centre %.0f kg, min over 12x8 m %.0f kg, area with >= 300 kg %.0f %%\n', ...
        zz(kz), pc, min(pz(:)), 100*mean(pz(:) >= 300));
end

%% 4. grid convergence (reference robot, z = 3 m)
Ng = [23 45 89 177];
pct = zeros(size(Ng));
for k = 1:length(Ng)
    nx = Ng(k); ny = round((Ng(k) - 1)*28/44) + 1;
    xg = linspace(-0.92*P.L, 0.92*P.L, nx);
    yg = linspace(-0.88*P.W, 0.88*P.W, ny);
    w = cdpr_wrench(P);
    c = 0;
    for i = 1:ny
        for j = 1:nx
            [~, ~, W] = cdpr_ik(P, [xg(j); yg(i); 3; 0; 0; 0]);
            c = c + cdpr_feasible(W, w, P.tmin, P.tmax);
        end
    end
    pct(k) = 100*c/(nx*ny);
    fprintf('grid %d x %d: %.2f %%\n', nx, ny, pct(k));
end

save results_validation.mat err_fk err_J agree Npt nfeas pay xs ys zz Ng pct

%% Figure: CoGiRo payload map at z = 1 m
figure('Position', [100 100 560 380]);
contourf(xs, ys, pay(:,:,1), 10); colorbar; axis equal tight;
hold on; plot(Pc.A(1,:), Pc.A(2,:), 'ks', 'MarkerFaceColor', 'k');
xlabel('x, m'); ylabel('y, m');
title('CoGiRo, z = 1 m: static payload capacity, kg');
print('-dpng', '-r200', 'fig_validation.png');
