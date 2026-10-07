function A = cdpr_tilt_all(P, G)
% pulley points when all 4 masts are tilted
% G - 4x2 matrix of tilt angles (rad): [outward in X, outward in Y] for each mast
% rotation is around the mast base (height P.zbase, default 0)

if ~isfield(P, 'zbase'), P.zbase = 0; end
A = P.A;
for j = 1:4
    sx = P.towers(j,1);
    sy = P.towers(j,2);
    c = [sx*P.L; sy*P.W; P.zbase];     % mast base
    a = G(j,1); b = G(j,2);
    Ry = [cos(a) 0 sx*sin(a); 0 1 0; -sx*sin(a) 0 cos(a)];
    Rx = [1 0 0; 0 cos(b) sy*sin(b); 0 -sy*sin(b) cos(b)];
    Rt = Ry*Rx;
    for i = [2*j-1, 2*j]
        A(:,i) = c + Rt*(P.A(:,i) - c);
    end
end
end
