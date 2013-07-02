function ErrXV = fsolve_ECI2Osc(Osc, X, V)
            % Osc = [e; p; theta; omega; Omega; i];
            
            [Xcalc, Vcalc] = Osc2ECI(Osc);
            
            ErrX = X - Xcalc;
            ErrV = V - Vcalc;
            
            ErrXV = [ErrX; ErrV];
            
end

