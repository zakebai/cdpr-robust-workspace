function s = cdpr_capacity(W, w, tmin, tmax)
% capacity margin by hyperplane shifting method
% (second method, does not use linear programming)
% s >= 0 : wrench w is inside the available wrench set
% s < 0  : outside

[m, n] = size(W);
C = nchoosek(1:n, m-1);       % groups of m-1 cables
s = inf;
for k = 1:size(C,1)
    N = null(W(:, C(k,:))');  % normal of the facet
    if size(N,2) ~= 1
        continue              % columns not independent
    end
    N = N/norm(N);
    p = N'*W;                 % projection of each cable column
    for sg = [1 -1]
        q = sg*p;
        tl = tmin(:)'.*ones(1,n);                 % tmin: one number or n numbers
        hmax = sum(q.*tl) + sum(max(q,0).*(tmax - tl));   % shifted hyperplane
        d = hmax - sg*(N'*w);
        s = min(s, d);
    end
end
end
