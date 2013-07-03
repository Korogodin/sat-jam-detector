function Osc = ECI2Osc(X, V, Osc0)
% Osc = [e; p; theta; omega; Omega; i];

% X = [x; y; z]
% V = [Vx; Vy; Vz]
options_solve = optimset('Display','off');  % Turn off display for fsolve
options_solve.MaxFunEvals = 6000;
options_solve.MaxIter = 4000;

for i = 1:4
    [Osc, fval, exitflag] = fsolve(@(OscVar)(fsolve_ECI2Osc(OscVar, X, V)), Osc0, options_solve);
    if exitflag ~= -2
        return;
    else
        Osc = Osc + [0; 0; 1; 1; 0; 0];
    end
end

fprintf('fsolve failed \n');

end

