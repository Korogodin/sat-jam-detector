function Osc = ECI2Osc(X, V, Osc0)
% Osc = [e; p; theta; omega; Omega; i];

% X = [x; y; z]
% V = [Vx; Vy; Vz]
options_solve = optimset('Display','off');  % Turn off display for fsolve
options_solve.MaxFunEvals = 6000;
options_solve.MaxIter = 4000;
mult_matrix = [1; 1; 1; 1; 1; 1];
Osc = fsolve(@(OscVar)(fsolve_ECI2Osc(OscVar./mult_matrix, X, V)), Osc0./mult_matrix, options_solve);
Osc = Osc .* mult_matrix;

end

