function [x, iter, res] = cdpr_fk(P, l_meas, x)
% forward kinematics, Gauss-Newton method with step halving
% (step is made smaller if the error does not decrease)
% l_meas - measured lengths, x - start pose

x = x(:);
h = 1e-7;
f = cdpr_ik(P, x) - l_meas(:);
for iter = 1:100
    % jacobian (numerical)
    J = zeros(8,6);
    for j = 1:6
        dx = zeros(6,1);
        dx(j) = h;
        J(:,j) = (cdpr_ik(P, x+dx) - cdpr_ik(P, x-dx)) / (2*h);
    end

    step = -(J'*J) \ (J'*f);
    a = 1;
    for k = 1:30
        fn = cdpr_ik(P, x + a*step) - l_meas(:);
        if norm(fn) < norm(f)
            break
        end
        a = a/2;
    end
    if norm(fn) >= norm(f)
        break                        % no better point: minimum found
    end
    x = x + a*step;
    f = fn;
    if norm(a*step) < 1e-10
        break
    end
end
res = norm(f);
end
