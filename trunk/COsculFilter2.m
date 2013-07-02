classdef COsculFilter2 < handle
    % X = [e; e'; p/p_mult; p'/p_mult; theta; theta'; omega; omega'; Omega; Omega'; i; i'];
    
    properties
        F
        T
        Dg
        Dn
        G
        c
        p_mult
        GDgG;   
        
        Xextr
        Xest
        Dest
        Dextr
        S
        
        Xforextr % [m]
        Vforextr % [m/s]
        Xforest % [m]
        Vforest % [m/s]       
        
        dY
        ImitJamAlert
        goodMeasLine;
        
        mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
        
        PredLe
        Xpred
        Ypred
        Zpred
    end
    
    methods
        
        function OF = COsculFilter2(T, Xest, Dest, p_mult, Dg, Dn, F, G, c)
            OF.goodMeasLine = 0;
            
            OF.T = T;
            
            OF.Xest = Xest;
            OF.Dest = Dest;
                        
            OF.p_mult = p_mult;
            
            OF.Dg = Dg;
            OF.GDgG =  G*Dg*G';   
            
            OF.Dn = Dn;
            
            OF.F = F;
            OF.G = G;
            OF.c = c;
            
            OF.PredLe = 20; 
            OF.Xpred = zeros(1, OF.PredLe);
            OF.Ypred = zeros(1, OF.PredLe);
            OF.Zpred = zeros(1, OF.PredLe);
        end
        
        function Extrapolate(OF)
            OF.Xextr = OF.F * OF.Xest;
            
           
%             OF.S = [1 0 0 0 0 0 0 0 0 0 0 0;
%                         0 0 0 0 0 0 0 0 0 0 0 0;
%                         0 0 1 0 0 0 0 0 0 0 0 0;
%                         0 0 0 0 0 0 0 0 0 0 0 0;
%                         0 0 0 0 1 0 0 0 0 0 0 0;
%                         0 0 0 0 0 0 0 0 0 0 0 0;
%                         0 0 0 0 0 0 1 0 0 0 0 0;
%                         0 0 0 0 0 0 0 0 0 0 0 0;
%                         0 0 0 0 0 0 0 0 1 0 0 0;
%                         0 0 0 0 0 0 0 0 0 0 0 0;
%                         0 0 0 0 0 0 0 0 0 0 1 0;
%                         0 0 0 0 0 0 0 0 0 0 0 0];
            OF.S = [1 0 0 0 0 0 0 0 0 0 0 0;
                        0 0 1 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 1 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 1 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 1 0 0 0;
                        0 0 0 0 0 0 0 0 0 0 1 0];
                    
            OF.Dextr = OF.F*OF.Dest*OF.F' + OF.GDgG;
        end
        
        
        function Estimate(OF, Osc)
            
            Osc(2) = Osc(2) / OF.p_mult;
            dY = Osc - [OF.Xextr(1); OF.Xextr(3); OF.Xextr(5); OF.Xextr(7); OF.Xextr(9); OF.Xextr(11)];            
            OF.dY = dY;
            
            t1 = OF.S'/OF.Dn*OF.S;
            t2 = inv(OF.Dextr);
            OF.Dest = inv(t1 + t2);
            
            OF.Xest = OF.Xextr + OF.Dest*OF.S'/OF.Dn*dY;
            OF.calcXVest();
        end
        
        function [X, V] = calcXV(OF, e, p, theta, omega, Omega, i0)
            
            if e > 0
                 theta =  OF.mod_pm_pi(theta);
                 omega =  OF.mod_pm_pi(omega);
            else
                 e = -e;
                 theta =  OF.mod_pm_pi(-theta);
                 omega =  OF.mod_pm_pi(-omega);
            end
            Omega =  OF.mod_pm_pi(Omega);
            i0 = OF.mod_pm_pi(i0);
            
            %Calc vectors x, y, z and Vx, Vy, Vz for orbital elements
            munapi = sqrt(OF.mu_earth / p);
            Vr = munapi*e*sin(theta);
            Vu = munapi*(1+e*cos(theta));
            r = p / (1+e*cos(theta));
            u = theta + omega;
            
            xyz = OF.U3(-Omega)*OF.U1(-i0)*OF.U3(-u)*[r; 0; 0];
            
            x = xyz(1);
            y = xyz(2);
            z = xyz(3);
            
            Vx = Vr.*x./r - Vu.*(sin(u).*cos(Omega) + cos(u).*sin(Omega).*cos(i0));
            Vy = Vr.*y./r - Vu.*(sin(u).*sin(Omega) - cos(u).*cos(Omega).*cos(i0));
            Vz = Vr.*z./r + Vu.*cos(u).*sin(i0);
            
            X = [x; y; z];
            V = [Vx; Vy; Vz];
        end
        
        function [X, V] = calcXVest(OF)
            %Calc vectors x, y, z and Vx, Vy, Vz for orbital elements estimations
            [X, V]  = OF.calcXV(OF.Xest(1), ...
                        OF.Xest(3) * OF.p_mult, ...
                        OF.Xest(5), ...
                        OF.Xest(7), ...
                        OF.Xest(9), ...
                        OF.Xest(11));
            OF.Xforest = X;
            OF.Vforest = V;
        end        
        
        function [X, V] = calcXVextr(OF)
            %Calc vectors x, y, z and Vx, Vy, Vz for orbital elements extrapolations
            [X, V]  = OF.calcXV(OF.Xextr(1), ...
                        OF.Xextr(3) * OF.p_mult, ...
                        OF.Xextr(5), ...
                        OF.Xextr(7), ...
                        OF.Xextr(9), ...
                        OF.Xextr(11));
            OF.Xforextr = X;
            OF.Vforextr = V;
        end      
        
        function Prediction(OF, Xmeas)
            if (norm(Xmeas - [OF.Xpred(1); OF.Ypred(1); OF.Zpred(1)]) > 100) && (OF.goodMeasLine > OF.PredLe)
                OF.ImitJamAlert = 1;
            else
                OF.ImitJamAlert = 0;                
            end
%                 OF.ImitJamAlert = norm(OF.Xforest - [OF.Xpred(1); OF.Ypred(1); OF.Zpred(1)]);            
            
            Xest = OF.Xest;
            for i = 1:OF.PredLe
                Xest = OF.F*Xest;
            end
            [X, V]  = OF.calcXV(Xest(1), ...
                Xest(3) * OF.p_mult, ...
                Xest(5), ...
                Xest(7), ...
                Xest(9), ...
                Xest(11));
            OF.Xpred(1) = X(1);
            OF.Xpred = circshift(OF.Xpred, [0; -1]);
            OF.Ypred(1) = X(2);
            OF.Ypred = circshift(OF.Ypred, [0; -1]);          
            OF.Zpred(1) = X(3);
            OF.Zpred = circshift(OF.Zpred, [0; -1]);            
        end
        
        function [X, V] = Osc2Orbit(OF, Osc, K)
            for k = 1:K
            end
        end
        
        function Osc = ECI2Oscul(OF, X, V, Osc0)
            % Osc = [theta; omega; Omega; i; e; p]
            % X = [x; y; z]
            % V = [Vx; Vy; Vz]
            options_solve = optimset('Display','on');  % Turn off display for fsolve
            options_solve.MaxFunEvals = 60000;
            options_solve.MaxIter = 40000;
            Osc = fsolve(@(Xfs)(OF.fsolve_Kepler(Xfs, X(1), X(2), X(3),...
                V(1), V(2), V(3))), Osc0, options_solve);
        end
        
        function X = Osc2X(OF, Osc, p_mult)
            % Osc = [theta; omega; Omega; i; e; p]
            % X = [e; e'; p/p_mult; p'/p_mult; theta; theta'; omega; omega'; Omega; Omega'; i; i'];
            X = zeros(12, 1);
            Xs = Osc;
            X(1) = Xs(5);
            X(3) = Xs(6) / p_mult;
            X(5) = Xs(1);
            X(7) = Xs(2);
            X(9) = Xs(3);
            X(11) = Xs(4);
        end
        
        function Osc = X2Osc(OF, X, p_mult)
            % Osc = [theta; omega; Omega; i; e; p]
            % X = [e; e'; p/p_mult; p'/p_mult; theta; theta'; omega; omega'; Omega; Omega'; i; i'];
            Osc = zeros(6,1);
            Osc(1) = X(5);
            Osc(2) = X(7);
            Osc(3) = X(9);
            Osc(4) = X(11);
            Osc(5) = X(1);
            Osc(6) = X(3)*p_mult;
        end        
        
        function FastInit(OF, X, V)
            % X = [x; y; z]
            % V = [Vx; Vy; Vz]
            Osc0 = OF.X2Osc(OF.Xest, OF.p_mult);
            Osc = OF.ECI2Oscul(X, V, Osc0);
            OF.Xest = OF.Osc2X(Osc, OF.p_mult);
        end
       
        function M = U1(OF, x)
            M = [1      0       0;
                0      cos(x)  sin(x);
                0      -sin(x) cos(x)];
        end
        
      
        function M = U3(OF, x)
            M = [cos(x)     sin(x)  0;
                -sin(x)    cos(x)  0;
                0          0       1];
        end
        
        function y = mod_pm_pi(OF, x )
            %MOD_PM_PI mod [-pi; pi];
            y = mod(x + pi, 2*pi) - pi;
        end
        
        function Ysolve = fsolve_Kepler(OF, Xs, Xizm, Yizm, Zizm, VXizm, VYizm, VZizm)
            
            %%%%Xsolve = r, u, Omega, i, d_r, d_u
            % Xsolve = theta, omega_p, Omega, i, e, p
            theta = Xs(1);
            omega_p = Xs(2);
            Omega = Xs(3);
            i = Xs(4);
            e = Xs(5);
            p = Xs(6);
            
            [X, V] = OF.calcXV( e, p, theta, omega_p, Omega, i);
            
            ErrX = Xizm - X(1);
            ErrY = Yizm - X(2);
            ErrZ = Zizm - X(3);
            
            ErrVx = VXizm - V(1);
            ErrVy = VYizm - V(2);
            ErrVz = VZizm - V(3);
            
            Ysolve = [ErrX, ErrY, ErrZ, ErrVx, ErrVy, ErrVz];
            
        end
        
    end
end

