function Wr = wind_poly(P, masses, Fw, n, outer)
% wrench vertices: mass interval x n-gon of wind forces
% outer = true: polygon around the disc |F| <= Fw (conservative)
% outer = false: polygon inside the disc (non-conservative)
if outer, r = Fw/cos(pi/n); else, r = Fw; end
ang = (0:n-1)*2*pi/n + pi/n;
Wr = [];
for m = masses
    for a = ang
        F = r*[cos(a); sin(a); 0];
        Wr = [Wr, [m*(-P.g(:)) - F; 0; 0; 0]];
    end
end
end
