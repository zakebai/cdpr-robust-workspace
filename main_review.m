%% Checks added after the supervisor review
% R1 computing time and number of LP problems
% R2 octagon vs 16-gon vs inscribed octagon (wind set approximation)
% R3 number of cable sample points in the clearance test (20 / 40 / 80)
% R4 pit robot: wind force doubled to include wind load on the cables
clear; clc; close all;
warning('off', 'all');
global NLP
NLP = 0;

P = cdpr_params(1, 'crossed');
rho = 1.225; Cd = 1.2; Ar = 0.96; v = 25;
Fw = 0.5*rho*Cd*Ar*v^2;

%% R1a: one LP test
rand('seed', 1);
N = 2000; t0 = tic;
for k = 1:N
    x = [12*(rand-0.5); 8*(rand-0.5); 0.5 + 5*rand; 0; 0; 0];
    [~, ~, W] = cdpr_ik(P, x);
    cdpr_feasible(W, cdpr_wrench(P), P.tmin, P.tmax, true);
end
t1 = toc(t0)/N;
fprintf('R1 one LP test (incl. IK): %.2f ms\n', 1000*t1);

%% tilt list (19 geometries), as in main_worstcase.m
g = deg2rad(1); Pt = {P};
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

nx = 60; dx = 2*P.L/nx; xs = -P.L + dx/2 : dx : P.L;
ny = round(nx*P.W/P.L); dy = 2*P.W/ny; ys = -P.W + dy/2 : dy : P.W;

%% R1b: full slice z = 3 m, level D (payload + wind + 19 tilts), counting LPs
WC = cdpr_wind_vertices(P, [350 650], Fw);
NLP = 0; t0 = tic; ok = 0;
for i = 1:length(ys)
    for j = 1:length(xs)
        ok = ok + robust_count(Pt, [xs(j); ys(i); 3], WC);
    end
end
tD = toc(t0);
fprintf('R1 slice z=3, level D: %.1f %% usable, %d points, %d LPs, %.1f s\n', ...
    100*ok/(length(xs)*length(ys)), length(xs)*length(ys), NLP, tD);

%% R2: wind polygon, level C (payload + wind, straight masts), z = 1..5
zs = 1:5;
Rw = zeros(length(zs), 3);
for k = 1:length(zs)
    for c = 1:3
        if c == 1, Wr = wind_poly(P, [350 650], Fw, 8, true);    % circumscribed octagon (paper)
        elseif c == 2, Wr = wind_poly(P, [350 650], Fw, 16, true); % circumscribed 16-gon
        else, Wr = wind_poly(P, [350 650], Fw, 8, false);        % inscribed octagon (non-conservative)
        end
        n = 0;
        for i = 1:length(ys)
            for j = 1:length(xs)
                n = n + cdpr_robust_ok({P}, [xs(j); ys(i); zs(k)], Wr);
            end
        end
        Rw(k,c) = 100*n/(length(xs)*length(ys));
    end
    fprintf('R2 z=%d: octagon %.1f | 16-gon %.1f | inscribed octagon %.1f %%\n', zs(k), Rw(k,:));
end

save results_review.mat t1 tD NLP Rw
end_of_ref = 1;

%% R3 and R4: pit robot, mast 50 m, sag 3 %
D = 100; cx = 140; cy = 90; alpha = deg2rad(55);
fx = cx - D/tan(alpha); fy = cy - D/tan(alpha);
zsurf = @(x, y) min(D, max(0, max(abs(x) - fx, abs(y) - fy))*tan(alpha));
P0 = cdpr_params(25, 'crossed'); P0.B = 0.2*P0.B; m0 = 50000;
P0.m = m0; P0.tmin = 0.02*m0*9.81; P0.tmax = 250e3; P0.rho_c = 0.0039*40^2; P0.zbase = D;
Hm = 50; P = P0; P.H = D + Hm; P.A(3,:) = D + Hm;
Ps = P; Ps.sag_eps = 0.03;
Pp = {Ps};
for j = 1:4
    for d = 1:2
        for sg = [1 -1]
            G = zeros(4,2); G(j,d) = sg*deg2rad(1);
            Q = Ps; Q.A = cdpr_tilt_all(Ps, G); Pp{end+1} = Q;
        end
    end
end
Q = Ps; Q.A = cdpr_tilt_all(Ps,  deg2rad(1)*ones(4,2)); Pp{end+1} = Q;
Q = Ps; Q.A = cdpr_tilt_all(Ps, -deg2rad(1)*ones(4,2)); Pp{end+1} = Q;
Fp = 0.5*1.225*1.2*24*25^2;
fprintf('cable wind per metre: %.1f N/m\n', 0.5*1.225*1.2*0.040*25^2);
h = 10; hc = 5;
xs = -cx + h/2 : h : cx; ys = -cy + h/2 : h : cy;
nfl = 0; rob = zeros(1,2); clr = zeros(1,3); npts = [20 40 80];
t0 = tic;
for i = 1:length(ys)
    for j = 1:length(xs)
        z0 = zsurf(xs(j), ys(i));
        if z0 > 0, continue; end            % floor only
        nfl = nfl + 1;
        x = [xs(j); ys(i); z0 + hc; 0; 0; 0];
        for c = 1:2
            WCp = cdpr_wind_vertices(P0, [0.7 1.3]*m0, c*Fp);
            rob(c) = rob(c) + cdpr_robust_ok(Pp, x(1:3), WCp);
        end
        [l] = cdpr_ik(P, x); Bw = x(1:3) + P.B;
        for q3 = 1:3
            ok = true;
            for c = 1:8
                lh = norm(P.A(1:2,c) - Bw(1:2,c)); fv = 0.03*l(c)^2/lh;
                for t = linspace(0.02, 0.98, npts(q3))
                    q = Bw(:,c) + t*(P.A(:,c) - Bw(:,c)); q(3) = q(3) - 4*fv*t*(1 - t);
                    if q(3) < zsurf(q(1), q(2)) + 1, ok = false; end
                end
            end
            clr(q3) = clr(q3) + ok;
        end
    end
end
fprintf('R4 pit floor robust (50 m, 3 %%): wind 11 kN %.0f %% | wind 22 kN %.0f %%\n', 100*rob/nfl);
fprintf('R3 floor cable clearance with 20/40/80 points: %.0f / %.0f / %.0f %%\n', 100*clr/nfl);
fprintf('pit part time %.0f s\n', toc(t0));
save -append results_review.mat rob clr nfl
