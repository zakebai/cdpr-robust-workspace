# Parabolic sag model (Step 7) against the exact inextensible catenary
# same end points, same horizontal tension H = tau*lh/l
# prints: sag error, error of the platform force (chord force + half cable weight), cable length - chord
import numpy as np
from scipy.optimize import brentq

q = 6.24*9.81          # 40 mm steel cable, weight per metre (N/m)

def catenary(lh, h, H):
    a = H/q
    x0 = brentq(lambda x0: a*(np.cosh((lh - x0)/a) - np.cosh(-x0/a)) - h, -50*lh, 50*lh)
    x = np.linspace(0, lh, 4001)
    y = a*np.cosh((x - x0)/a); y -= y[0]
    l = np.hypot(lh, h)
    sag = (x*h/l - y*lh/l).max()                       # largest distance from the chord
    arc = a*(np.sinh((lh - x0)/a) - np.sinh(-x0/a))    # cable length
    F = H*np.array([1, np.sinh(-x0/a)])                # force on the platform (horizontal, vertical)
    return sag, arc, F, l

print('lh   h  eps | sag err %  force err (% of tau)  length-chord: parabola / catenary (m)')
for lh, h in [(150, 120), (220, 145), (300, 150)]:
    for eps in [0.01, 0.03, 0.05, 0.08, 0.11]:
        l = np.hypot(lh, h)
        tau = q*lh/(8*eps)                             # tension that gives sag ratio eps (Eq. 8)
        fp = q*l*lh/(8*tau)                            # parabolic sag across the chord
        fc, arc, F, _ = catenary(lh, h, tau*lh/l)
        Fm = tau*np.array([lh/l, h/l]) - np.array([0, q*l/2])   # model force
        err = np.linalg.norm(F - Fm)/tau*100
        print(f'{lh} {h} {eps:.2f} | {100*(fc/fp - 1):+.1f}  {err:.2f}  {8*fp**2/(3*l):.2f} / {arc - l:.2f}')
