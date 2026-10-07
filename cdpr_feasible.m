function ok = cdpr_feasible(W, w, tmin, tmax, lp_only)
% check if tensions exist: W*tau = w, tmin <= tau <= tmax
% (linear programming, only feasibility)
% the found tensions are checked again; if the solver is not sure,
% the second method (capacity margin) decides
% lp_only = 1 : use only linear programming (for validation)

if nargin < 5, lp_only = 0; end
n = size(W,2);
% scale tensions to 0..1 (better numerics)
sc = tmax;
ws = w/sc;
lb = tmin(:)/sc.*ones(n,1);    % tmin can be one number or 8 numbers
ub = ones(n,1);
f = zeros(n,1);

if exist('linprog')
    opt = optimoptions('linprog', 'Display', 'off');
    [t, ~, flag] = linprog(f, [], [], W, ws, lb, ub, opt);
    solved = (flag == 1 || flag == -2);       % -2: no feasible point
    ok = (flag == 1);
else
    % Octave: glpk
    ctype = repmat('S', 1, size(W,1));
    vtype = repmat('C', 1, n);
    prm.msglev = 0;
    prm.itlim = 5000;          % rare degenerate problems can cycle
    [t, ~, err, extra] = glpk(f, W, ws, lb, ub, ctype, vtype, 1, prm);
    solved = (err == 0) || (err >= 10);       % >= 10: no feasible point
    ok = (err == 0) && (extra.status == 5 || extra.status == 2);
end

% check the found tensions (solver tolerance)
if ok
    good = norm(W*t - ws) < 1e-6 && all(t >= lb - 1e-6) && all(t <= ub + 1e-6);
    if ~good
        ok = false;
        solved = false;        % not sure
    end
end

if ~solved && ~lp_only
    ok = cdpr_capacity(W, w, tmin, tmax) >= 0;
end
end
