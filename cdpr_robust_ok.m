function ok = cdpr_robust_ok(Plist, x, Wr, ori)
% robust check: pose x is feasible for ALL cases
% Plist - cell array of robots (for example nominal and tilted masts)
% Wr    - 6xK external wrenches (vertices of mass and wind sets)
% ori   - 3xQ orientations [phi; theta; psi] to test at this point
% stops at the first failed case

if nargin < 4, ori = [0;0;0]; end
ok = true;
for q = 1:size(ori,2)
    xq = [x(1:3); ori(:,q)];
    for p = 1:numel(Plist)
        P = Plist{p};
        [l, ~, W] = cdpr_ik(P, xq);
        tmin = P.tmin;
        dw = zeros(6,1);
        if isfield(P, 'sag_eps')
            tmin = cdpr_tmin_sag(P, xq);    % sag limit for long cables
            dw(3) = P.rho_c*9.81*sum(l)/2;  % half of cable weight hangs on platform
        end
        for k = 1:size(Wr,2)
            if ~cdpr_feasible(W, Wr(:,k) + dw, tmin, P.tmax)
                ok = false;
                return
            end
        end
    end
end
end
