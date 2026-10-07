function tmin = cdpr_tmin_sag(P, x)
% minimum tension of each cable so that the sag is not more than
% P.sag_eps * cable length (parabolic cable, P.rho_c kg/m)
% sag across the chord at mid-span: f = rho*g*l*lh/(8*tau)
% (vertical sag is f*l/lh)
% f <= eps*l   ->   tau >= rho*g*lh/(8*eps)

R = rotm_zyx(x(4), x(5), x(6));
Bw = R*P.B + repmat(x(1:3), 1, 8);
Lv = P.A - Bw;
lh = sqrt(sum(Lv(1:2,:).^2, 1))';
tau_sag = P.rho_c*9.81*lh/(8*P.sag_eps);
tmin = max(P.tmin, tau_sag);
end
