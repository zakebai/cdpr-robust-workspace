function Wr = cdpr_wind_vertices(P, masses, Fw)
% wrench vertices for mass interval and horizontal wind |F| <= Fw
% the wind disc is covered by an octagon (radius Fw/cos(pi/8)),
% so checking the vertices is enough (feasible wrench set is convex)

r = Fw/cos(pi/8);
ang = (0:7)*pi/4 + pi/8;
Wr = [];
for m = masses
    for a = ang
        F = r*[cos(a); sin(a); 0];
        Wr = [Wr, [m*(-P.g(:)) - F; 0; 0; 0]];
    end
end
if Fw == 0
    Wr = [];
    for m = masses
        Wr = [Wr, [m*(-P.g(:)); 0; 0; 0]];
    end
end
end
