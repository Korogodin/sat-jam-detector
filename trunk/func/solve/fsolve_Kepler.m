function Ysolve = fsolve_Kepler(Xs, Xizm, Yizm, Zizm, VXizm, VYizm, VZizm)

global mu_earth

%%%%Xsolve = r, u, Omega, i, d_r, d_u
% Xsolve = theta, omega_p, Omega, i, e, p
theta = Xs(1);
omega_p = Xs(2);
Omega = Xs(3);
i = Xs(4);
e = Xs(5);
p = Xs(6);

[x, y, z, Vx, Vy, Vz] = get_vector_XV( e, p, theta, omega_p, Omega, i);

ErrX = Xizm - x;
ErrY = Yizm - y;
ErrZ = Zizm - z;

ErrVx = VXizm - Vx;
ErrVy = VYizm - Vy;
ErrVz = VZizm - Vz;

Ysolve = [ErrX, ErrY, ErrZ, ErrVx, ErrVy, ErrVz];

end

