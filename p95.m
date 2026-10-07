function q = p95(v)
% 95th percentile (simple, no toolbox)
v = sort(v(:));
q = v(max(1, ceil(0.95*length(v))));
end
