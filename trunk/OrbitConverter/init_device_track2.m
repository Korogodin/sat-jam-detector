% Filter step
T = 1;
dTmod = T;

addpath('../');
const_OrbitConverter;

% Evolutionary matrix
F = [1   T   0   0   0   0   0   0   0   0   0   0;
       0   1   0   0   0   0   0   0   0   0   0   0;
       0   0   1   T   0   0   0   0   0   0   0   0;
       0   0   0   1   0   0   0   0   0   0   0   0;
       0   0   0   0   1   T   0   0   0   0   0   0;
       0   0   0   0   0   1   0   0   0   0   0   0;
       0   0   0   0   0   0   1   T   0   0   0   0;
       0   0   0   0   0   0   0   1   0   0   0   0;
       0   0   0   0   0   0   0   0   1   T   0   0;
       0   0   0   0   0   0   0   0   0   1   0   0;
       0   0   0   0   0   0   0   0   0   0   1   T;
       0   0   0   0   0   0   0   0   0   0   0   1];

p_mult = 5e7; % To simplify matrix calculations - reducing the dynamic range

Osc = [0; 7e6; 0; 0; 0; 0.85];
Xest =  [Osc(1); 0; Osc(2) / p_mult; 0; Osc(3); 0; Osc(4); 0; Osc(5); 0; Osc(6); 0];
Xextr = Xest;

% RMS of shaping noises
std_e = 5e-7 / 15*dTmod;
std_p = 5 / 15*dTmod / p_mult; % [m]
std_theta = 1e-5 / 15*dTmod; % [rad]
std_omega = 1e-9 / 15*dTmod; % [rad]
std_Omega = 1e-7 / 15*dTmod; % [rad]
std_i = 1e-8 / 15*dTmod; % [rad]

Dest = [std_e^2*1e1     0           0           0               0               0               0               0               0               0               0           0
            0           std_e^2*1e3 0           0               0               0               0               0               0               0               0           0
            0           0           std_p^2*1e3 0               0               0               0               0               0               0               0           0
            0           0           0           std_p^2*1e5     0               0               0               0               0               0               0           0
            0           0           0           0               std_theta^2*1e2 0               0               0               0               0               0           0
            0           0           0           0               0               std_theta^2*1e3 0               0               0               0               0           0
            0           0           0           0               0               0               std_omega^2*1e2 0               0               0               0           0
            0           0           0           0               0               0               0               std_omega^2*1e3 0               0               0           0
            0           0           0           0               0               0               0               0               std_Omega^2*1e2 0               0           0
            0           0           0           0               0               0               0               0               0               std_Omega^2*1e3 0           0
            0           0           0           0               0               0               0               0               0               0               std_i^2*1e2 0
            0           0           0           0               0               0               0               0               0               0               0           std_i^2*1e3];

G =  [0 0 0 0 0 0;
        1 0 0 0 0 0;
        0 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 0;
        0 0 0 0 0 1];

Dg = [std_e^2 0       0           0           0           0
      0       std_p^2 0           0           0           0 
      0       0       std_theta^2 0           0           0
      0       0       0           std_omega^2 0           0
      0       0       0           0           std_Omega^2 0
      0       0       0           0           0           std_i^2];

GDgG = G*Dg*G';                        

Dn = [1e-8   0      0    0       0       0
         0   1e-10      0    0       0       0
         0   0      1e-6    0       0       0
         0   0      0    1e-14       0       0
         0   0      0    0       1e-10       0
         0   0      0    0       0       1e-12];
 
c = [1 0 0 0 0 0 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 1 0 0 0 0 0 0 0;
      0 0 0 0 0 0 1 0 0 0 0 0;
      0 0 0 0 0 0 0 0 1 0 0 0; 
      0 0 0 0 0 0 0 0 0 0 1 0];
 
OF = COsculFilter2(T, Xest, Dest, p_mult, Dg, Dn, F, G, c);

