function ok = robust_count(Plist, x, Wr)
% same as cdpr_robust_ok (zero orientation, no sag), but counts LP problems
global NLP
ok = true;
for p = 1:numel(Plist)
    P = Plist{p};
    [~, ~, W] = cdpr_ik(P, [x(:); 0; 0; 0]);
    for k = 1:size(Wr,2)
        NLP = NLP + 1;
        if ~cdpr_feasible(W, Wr(:,k), P.tmin, P.tmax)
            ok = false; return
        end
    end
end
end
