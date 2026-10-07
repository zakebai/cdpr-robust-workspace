function mmax = cdpr_maxload(P, x, mmax_try)
% largest total mass (kg) the robot can hold at pose x
% (LP: maximize m, W*tau = m*wg, tmin <= tau <= tmax)
% the centre of mass P.com is given in platform frame

if nargin < 3, mmax_try = 1e5; end
if ~isfield(P, 'com'), P.com = [0;0;0]; end

[~, ~, W] = cdpr_ik(P, x);
R = rotm_zyx(x(4), x(5), x(6));
g = P.g(:);
c = R*P.com(:);
wg = [-g; -cross(c, g)];     % wrench for 1 kg

sc = P.tmax;
n = 8;
ms = sc/9.81;                % mass scale
Aeq = [W, -wg*ms/sc];        % variables: tau/sc (8) and m/ms
beq = zeros(6,1);
lb = [P.tmin/sc*ones(n,1); 0];
ub = [ones(n,1); mmax_try/ms];
f = [zeros(n,1); -1];

if exist('linprog')
    opt = optimoptions('linprog', 'Display', 'off');
    [v, ~, flag] = linprog(f, [], [], Aeq, beq, lb, ub, opt);
    ok = (flag == 1);
else
    ctype = repmat('S', 1, 6);
    vtype = repmat('C', 1, n+1);
    prm.msglev = 0;
    [v, ~, err, extra] = glpk(f, Aeq, beq, lb, ub, ctype, vtype, 1, prm);
    ok = (err == 0) && (extra.status == 5 || extra.status == 2);
end

if ok
    mmax = v(end)*ms;
else
    mmax = 0;
end
end
