function e = cdpr_eps_min(P, x, w)
% smallest sag ratio (sag / cable length) that the suspended robot can keep
% at pose x with wrench w. Bisection: smaller sag needs larger tmin,
% but in a suspended robot all tensions are limited by the weight.

lo = 1e-3; hi = 0.5;
[l, ~, W] = cdpr_ik(P, x);
dw = zeros(6,1); dw(3) = P.rho_c*9.81*sum(l)/2;   % half of cable weight
P.sag_eps = hi;
if ~cdpr_feasible(W, w + dw, cdpr_tmin_sag(P, x), P.tmax)
    e = inf;
    return
end
for k = 1:30
    mid = sqrt(lo*hi);
    P.sag_eps = mid;
    if cdpr_feasible(W, w + dw, cdpr_tmin_sag(P, x), P.tmax)
        hi = mid;
    else
        lo = mid;
    end
end
e = hi;
end
