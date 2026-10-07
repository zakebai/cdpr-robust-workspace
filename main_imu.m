%% Mast tilt: positioning error without and with IMU correction
% all 4 masts tilt randomly up to 1 deg (X and Y)
% the robot measures cable lengths and finds its pose by forward kinematics
% without correction: model with straight masts
% with correction: model uses tilt measured by IMU (error +- delta)

clear; clc; close all;
warning('off', 'all');
rng(2);

P = cdpr_params(1, 'crossed');
gmax = deg2rad(1);
delta = deg2rad([0.01 0.05 0.1]);
Ns = 200;          % tilt samples
Np = 5;            % poses per sample

e0 = zeros(Ns*Np, 1);
e1 = zeros(Ns*Np, length(delta));
o0 = zeros(Ns*Np, 1);
o1 = zeros(Ns*Np, length(delta));
n = 0;
for s = 1:Ns
    G = gmax*(2*rand(4,2) - 1);
    Pt = P; Pt.A = cdpr_tilt_all(P, G);         % real robot
    for p = 1:Np
        x = [-4 + 8*rand; -2.5 + 5*rand; 1 + 3*rand; 0; 0; 0];
        l = cdpr_ik(Pt, x);                     % measured lengths
        n = n + 1;
        xe = cdpr_fk(P, l, x);                  % no correction
        e0(n) = 1000*norm(xe(1:3) - x(1:3));
        o0(n) = rad2deg(max(abs(mod(xe(4:6) + pi, 2*pi) - pi)));
        for d = 1:length(delta)
            Gm = G + delta(d)*(2*rand(4,2) - 1);    % IMU reading
            Pm = P; Pm.A = cdpr_tilt_all(P, Gm);
            xe = cdpr_fk(Pm, l, x);
            e1(n,d) = 1000*norm(xe(1:3) - x(1:3));
            o1(n,d) = rad2deg(max(abs(mod(xe(4:6) + pi, 2*pi) - pi)));
        end
    end
end

fprintf('No correction: mean %.0f mm, 95%% %.0f mm, max %.0f mm | orientation 95%% %.2f deg, max %.2f deg\n', ...
    mean(e0), p95(e0), max(e0), p95(o0), max(o0));
for d = 1:length(delta)
    fprintf('IMU +-%.2f deg: mean %.1f mm, 95%% %.1f mm, max %.1f mm | orientation 95%% %.2f deg\n', ...
        rad2deg(delta(d)), mean(e1(:,d)), p95(e1(:,d)), max(e1(:,d)), p95(o1(:,d)));
end
save results_imu.mat e0 e1 o0 o1 delta

figure('Position', [100 100 500 350]);
bar([p95(e0), p95(e1(:,1)), p95(e1(:,2)), p95(e1(:,3))]); set(gca, 'YScale', 'log');
set(gca, 'XTickLabel', {'no IMU', '0.01', '0.05', '0.1'});
ylabel('position error, mm');
print('-dpng', '-r200', 'fig_imu.png');
