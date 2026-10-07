%% Application to a representative small deep open pit
% pit: crest 280 x 180 m, depth 100 m, overall slope 55 deg
% robot: 4 masts on the rim (300 x 200 m), mast height Hm above the rim
% part 1: smallest possible cable sag (suspended robot, tension limited by weight)
% part 2: checks on the pit surface:
%   1 rigid  - straight cables, nominal load (classic model)
%   2 robust - mass +-30 %, wind 25 m/s, mast tilt 1 deg, sag <= eps of
%              cable length, half of cable weight on the platform
%   3 robust + clearance of cables (with sag) and platform from the ground

clear; clc; close all;
warning('off', 'all');

%% Pit
D = 100;                             % depth, m
cx = 140; cy = 90;                   % crest half sizes, m
alpha = deg2rad(55);                 % overall slope
fx = cx - D/tan(alpha);              % floor half sizes
fy = cy - D/tan(alpha);
zsurf = @(x, y) min(D, max(0, max(abs(x) - fx, abs(y) - fy))*tan(alpha));
fprintf('pit floor %.0f x %.0f m\n', 2*fx, 2*fy);

%% Robot (industrial size)
P0 = cdpr_params(25, 'crossed');     % masts at x = +-150, y = +-100 m
P0.B = 0.2*P0.B;                     % platform 8 x 5 x 3 m
m0 = 50000;                          % kg, suspended mass (platform + payload)
P0.m = m0;
P0.tmin = 0.02*m0*9.81;
P0.tmax = 250e3;                     % N, 40 mm steel rope, safety factor ~4
P0.rho_c = 0.0039*40^2;              % kg/m
P0.zbase = D;                        % masts stand on the rim

%% Part 1: smallest sag ratio at the floor centre and floor corner
mlist = [10 20 35 50 80]*1000;
Hm_list = [50 100 150];
EPS = zeros(length(Hm_list), length(mlist), 2);
for kh = 1:length(Hm_list)
    P = P0; P.A(3,:) = D + Hm_list(kh);
    for km = 1:length(mlist)
        P.m = mlist(km); P.tmin = 0.02*P.m*9.81;
        EPS(kh,km,1) = cdpr_eps_min(P, [0; 0; 5; 0; 0; 0], cdpr_wrench(P));
        EPS(kh,km,2) = cdpr_eps_min(P, [fx; fy; 5; 0; 0; 0], cdpr_wrench(P));
        % closed form, Eq. (9), at the floor centre
        [l, u] = cdpr_ik(P, [0; 0; 5; 0; 0; 0]);
        Lv = P.A - ([0; 0; 5] + P.B);
        lh = sqrt(sum(Lv(1:2,:).^2, 1))';
        meff = P.m + P.rho_c*sum(l)/2;
        EPS(kh,km,3) = P.rho_c*sum(lh.*u(3,:)')/(8*meff);
    end
end
disp('smallest sag ratio, % (rows: mast 50/100/150 m; columns: 10/20/35/50/80 t)');
disp('floor centre'); disp(round(1000*EPS(:,:,1))/10);
disp('floor corner'); disp(round(1000*EPS(:,:,2))/10);
disp('floor centre, closed form Eq. (9)'); disp(round(1000*EPS(:,:,3))/10);

Ar = 24; Fw = 0.5*1.225*1.2*Ar*25^2;
WC = cdpr_wind_vertices(P0, [0.7 1.3]*m0, Fw);
fprintf('wind force %.1f kN\n', Fw/1000);

h = 10;                              % grid step, m
hc = 5;                              % platform centre above ground, m
xs = -cx + h/2 : h : cx;
ys = -cy + h/2 : h : cy;
ZZ = zeros(length(ys), length(xs));
for i = 1:length(ys)
    for j = 1:length(xs)
        ZZ(i,j) = zsurf(xs(j), ys(i));
    end
end
floor_ = (ZZ == 0);
wall = (ZZ > 0) & (ZZ < D);
lev = 0:20:80;

%% Part 2: pit surface for design cases [mast height, allowed sag ratio]
CASES = [50 0.03; 50 0.05; 100 0.03; 100 0.05];
RES = zeros(size(CASES,1), 8);
COV = zeros(size(CASES,1), length(lev));
MAPS = cell(size(CASES,1), 1);
for kh = 1:size(CASES,1)
    Hm = CASES(kh,1);
    eps_sag = CASES(kh,2);
    P = P0;
    P.H = D + Hm;
    P.A(3,:) = D + Hm;
    Ps = P; Ps.sag_eps = eps_sag;

    % robust set: straight masts + each mast tilted 1 deg in +-X, +-Y
    Pt = {Ps};
    for j = 1:4
        for d = 1:2
            for sg = [1 -1]
                G = zeros(4,2); G(j,d) = sg*deg2rad(1);
                Q = Ps; Q.A = cdpr_tilt_all(Ps, G); Pt{end+1} = Q;
            end
        end
    end
    Q = Ps; Q.A = cdpr_tilt_all(Ps,  deg2rad(1)*ones(4,2)); Pt{end+1} = Q;   % all out
    Q = Ps; Q.A = cdpr_tilt_all(Ps, -deg2rad(1)*ones(4,2)); Pt{end+1} = Q;   % all in

    M = zeros(size(ZZ));             % 0 out, 1 rigid, 2 robust, 3 robust + clearance
    for i = 1:length(ys)
        for j = 1:length(xs)
            x = [xs(j); ys(i); ZZ(i,j) + hc; 0; 0; 0];
            [l, ~, W] = cdpr_ik(P, x);
            if ~cdpr_feasible(W, cdpr_wrench(P), P.tmin, P.tmax), continue; end
            M(i,j) = 1;
            if ~cdpr_robust_ok(Pt, x(1:3), WC), continue; end
            M(i,j) = 2;
            ok = true;
            Bw = x(1:3) + P.B;
            for c = 1:8
                lh = norm(P.A(1:2,c) - Bw(1:2,c));
                fv = eps_sag*l(c)*l(c)/lh;          % vertical sag (eps is sag across the chord)
                for t = linspace(0.02, 0.98, 40)
                    q = Bw(:,c) + t*(P.A(:,c) - Bw(:,c));
                    q(3) = q(3) - 4*fv*t*(1 - t);     % sag shape
                    if q(3) < zsurf(q(1), q(2)) + 1, ok = false; end
                end
            end
            for sx = [-1 1]
                for sy = [-1 1]
                    q = x(1:3) + [sx*4; sy*2.5; -1.5];
                    if q(3) < zsurf(q(1), q(2)) + 1, ok = false; end
                end
            end
            if ok, M(i,j) = 3; end
        end
    end
    MAPS{kh} = M;
    RES(kh,:) = [Hm, 100*eps_sag, 100*mean(M(floor_) >= 1), 100*mean(M(floor_) >= 2), 100*mean(M(floor_) >= 3), ...
        100*mean(M(wall) >= 1), 100*mean(M(wall) >= 2), 100*mean(M(wall) >= 3)];
    for k = 1:length(lev)
        if k == 1
            sel = floor_;
        else
            sel = (ZZ >= lev(k-1) + 1e-9) & (ZZ < lev(k));
        end
        COV(kh,k) = 100*mean(M(sel) >= 3);
    end
    fprintf('Hm = %3d m, sag %.0f %% | floor: rigid %3.0f, robust %3.0f, +clearance %3.0f %% | walls: %3.0f, %3.0f, %3.0f %%\n', RES(kh,:));
end

disp('usable % by level (floor, 0-20, 20-40, 40-60, 60-80 m above floor)');
disp([CASES COV]);

P = P0; P.A(3,:) = D + 50;
[l, ~, ~] = cdpr_ik(P, [0; 0; hc; 0; 0; 0]);
fprintf('floor centre, mast 50 m: cable length %.0f m, cable weight on platform %.1f t\n', ...
    mean(l), P.rho_c*sum(l)/2/1000);

save results_pit.mat xs ys ZZ MAPS RES COV lev CASES EPS mlist Hm_list D fx fy cx cy Fw

figure('Position', [100 100 1000 600]);
for kh = 1:size(CASES,1)
    subplot(2, 2, kh);
    imagesc(xs, ys, MAPS{kh}); axis equal tight; set(gca, 'YDir', 'normal'); caxis([0 3]);
    hold on; rectangle('Position', [-fx -fy 2*fx 2*fy], 'EdgeColor', 'w');
    title(sprintf('mast %d m, sag %.0f %%', CASES(kh,1), 100*CASES(kh,2)));
end
print('-dpng', '-r200', 'fig_pit.png');
