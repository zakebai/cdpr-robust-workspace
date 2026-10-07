function P = cogiro_params()
% CoGiRo robot (LIRMM / Tecnalia), measured geometry
% source: CoGiRo description, Table I (November 2015), CASPR model files
% frame 15.24 x 11.24 x 5.93 m, 8 cables, tmax = 5000 N

P.A = [-7.1775 -5.4361 5.3911;
       -7.4594 -5.1504 5.3999;
       -7.3911  5.1940 5.3976;
       -7.1026  5.4753 5.4094;
        7.2398  5.3759 5.4093;
        7.5208  5.0851 5.4200;
        7.4461 -5.2539 5.3874;
        7.1608 -5.5342 5.3973]';

P.B = [ 0.5032 -0.4928  0.0;
       -0.5097  0.3508  0.9976;
       -0.5032 -0.2700  0.0;
        0.4960  0.3561  0.9996;
       -0.5032  0.4928  0.0;
        0.4998 -0.3404  0.9991;
        0.5021  0.2750 -0.0007;
       -0.5045 -0.3463  0.9976]';

P.m = 91.058;                      % platform only
P.com = [-0.034; -0.013; 0.264];   % centre of mass, platform frame
P.g = [0; 0; -9.81];
P.tmin = 100;
P.tmax = 5000;
end
