function [X, V] = Osc2ECI( Osc )
% Osc = [e; p; theta; omega; Omega; i];

    const_OrbitConverter;

    e = Osc(1);
    p = Osc(2);
    theta = Osc(3);
    omega = Osc(4);
    Omega = Osc(5);
    i0 = Osc(6);

    if e > 0
        theta =  mod_pm_pi(theta);
        omega = mod_pm_pi(omega);
    else
        e = -e;
        theta =  mod_pm_pi(-theta);
        omega =  mod_pm_pi(-omega);
    end
    Omega =  mod_pm_pi(Omega);
    i0 = mod_pm_pi(i0);

    %Calc vectors x, y, z and Vx, Vy, Vz for orbital elements
    munapi = sqrt(mu_earth / p);
    Vr = munapi*e*sin(theta);
    Vu = munapi*(1+e*cos(theta));
    r = p / (1+e*cos(theta));
    u = theta + omega;

    xyz = U3(-Omega)*U1(-i0)*U3(-u)*[r; 0; 0];

    x = xyz(1);
    y = xyz(2);
    z = xyz(3);

    Vx = Vr.*x./r - Vu.*(sin(u).*cos(Omega) + cos(u).*sin(Omega).*cos(i0));
    Vy = Vr.*y./r - Vu.*(sin(u).*sin(Omega) - cos(u).*cos(Omega).*cos(i0));
    Vz = Vr.*z./r + Vu.*cos(u).*sin(i0);

    X = [x; y; z];
    V = [Vx; Vy; Vz];


end

