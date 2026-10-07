# Robust workspace calculation for outdoor cable-driven robots

MATLAB / GNU Octave code for the article
*"Robust workspace calculation for outdoor cable-driven robots considering structural and environmental disturbances"* (MethodsX, submitted).

The code computes the wrench-feasible workspace of an eight-cable suspended cable-driven parallel robot (CDPR) when payload variation, wind from any horizontal direction, mast (anchor) tilt, cable sag and terrain clearance act together.

## Requirements
- MATLAB R2020b or later with Optimization Toolbox (`linprog`), **or**
- GNU Octave 8 or later (uses the built-in `glpk`).
- `catenary_check.py` needs Python 3 with NumPy and SciPy (optional).

## Quick start
```matlab
P  = cdpr_params(1, 'crossed');               % reference robot (12 x 8 x 6 m)
Fw = 0.5*1.225*1.2*0.96*25^2;                 % wind force bound, N
Wr = cdpr_wind_vertices(P, [350 650], Fw);    % 6 x 16 wrench vertices
G  = zeros(4,2); G(1,1) = deg2rad(1);         % mast 1 tilted 1 deg outward
Q  = P; Q.A = cdpr_tilt_all(P, G);            % tilted geometry
ok = cdpr_robust_ok({P, Q}, [0; 0; 3], Wr)    % robust test at one point
```

## Functions
| Function | Purpose | Eq. in article |
|---|---|---|
| `cdpr_params`, `cogiro_params` | Robot parameters (reference robot, measured CoGiRo geometry) | Table 2 |
| `cdpr_ik`, `cdpr_fk`, `rotm_zyx` | Inverse/forward kinematics, structure matrix W | (1)–(2) |
| `cdpr_wrench` | External wrench (gravity, wind) | (3) |
| `cdpr_feasible` | LP wrench-feasibility test | (4) |
| `cdpr_capacity` | Hyperplane-shifting capacity margin | (5) |
| `cdpr_wind_vertices` | Payload x wind octagon wrench vertices | (6) |
| `cdpr_tilt_all` | Anchor points of tilted masts | (7) |
| `cdpr_tmin_sag`, `cdpr_eps_min` | Sag-dependent minimum tension; smallest sag ratio | (8)–(9) |
| `cdpr_robust_ok` | Robust test over geometries, wrenches, orientations | (4), (6)–(8) |
| `cdpr_maxload` | Largest static payload at a pose | (4) |

## Scripts that reproduce the article
| Script | Results |
|---|---|
| `main_validation.m` | Validation V1–V4 (Table 4, Fig. 3) |
| `main_workspace3d.m` | Case study 1 (Table 5, Fig. 4) |
| `main_worstcase.m` | Case study 2 (Table 6, Fig. 5) |
| `main_imu.m` | Case study 3 (Fig. 6) |
| `main_pit.m` | Case study 4 (Tables 7–8, Fig. 7) |
| `main_review.m`, `main_review_tilt.m`, `catenary_check.py` | Computing cost, wind polygon, clearance sampling, cable wind, random tilts, catenary check |

Folder `old/` contains earlier scripts that are not used in the article.

## License
MIT (see `LICENSE`).

## How to cite
See `CITATION.cff`.
