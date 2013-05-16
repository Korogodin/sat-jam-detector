classdef COsculFilter < handle
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
        
        mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
        
        PredLe
        Xpred
        Ypred
        Zpred
    end
    
    methods
        
        function OF = COsculFilter(T, Xest, Dest, p_mult, Dg, Dn, F, G, c)
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
            
            e = OF.Xextr(1);
            p = OF.Xextr(3) * OF.p_mult;
            theta = OF.Xextr(5);
            omega = OF.Xextr(7);
            Omega = OF.Xextr(9);
            i0 = OF.Xextr(11);

            OF.calcXVextr;

            munapi = sqrt(OF.mu_earth / p);
            u = theta + omega;
            
            oec = 1 + e*cos(theta);
            es = e*sin(theta);
            
            Sc_x = cos(u)*cos(Omega) - sin(u)*sin(Omega)*cos(i0);
            Sc_y = cos(u)*sin(Omega) + sin(u)*cos(Omega)*cos(i0);
            Sc_z = sin(u)*sin(i0);
            
            Sc_Vx = sin(u)*cos(Omega) + cos(u)*sin(Omega)*cos(i0);
            Sc_Vy = sin(u)*sin(Omega) - cos(u)*cos(Omega)*cos(i0);
            Sc_Vz = cos(u)*sin(i0);
            
            S0 = nan(6,6);
            
            % Discriminator for e
            S0(1, 1) = - p * Sc_x / (1+e*cos(theta))^2 * cos(theta);
            S0(2, 1) = - p * Sc_y / (1+e*cos(theta))^2 * cos(theta);
            S0(3, 1) = - p * Sc_z / (1+e*cos(theta))^2 * cos(theta);
            S0(4, 1) = munapi * (sin(theta)*Sc_x - cos(theta)*Sc_Vx);
            S0(5, 1) = munapi * (sin(theta)*Sc_y - cos(theta)*Sc_Vy);
            S0(6, 1) = munapi * (sin(theta)*Sc_z + cos(theta)*Sc_Vz);
            
            % Discriminator for p
            S0(1, 2) = Sc_x/oec;
            S0(2, 2) = Sc_y/oec;
            S0(3, 2) = Sc_z/oec;
            S0(4, 2) = - 0.5 * 1/p * munapi * (es *Sc_x - oec*Sc_Vx);
            S0(5, 2) = - 0.5 * 1/p * munapi * (es *Sc_y - oec*Sc_Vy);
            S0(6, 2) = - 0.5 * 1/p * munapi * (es *Sc_z + oec*Sc_Vz);
            S0(:,2) = S0(:,2) * OF.p_mult;
            
            % Discriminator for theta
            Sc_x_theta = - sin(u)*cos(Omega) - cos(u)*sin(Omega)*cos(i0);
            Sc_y_theta = - sin(u)*sin(Omega) + cos(u)*cos(Omega)*cos(i0);
            Sc_z_theta = cos(u)*sin(i0);
            Sc_Vx_theta = cos(u)*cos(Omega) - sin(u)*sin(Omega)*cos(i0);
            Sc_Vy_theta = cos(u)*sin(Omega) + sin(u)*cos(Omega)*cos(i0);
            Sc_Vz_theta = -sin(u)*sin(i0);
            
            S0(1, 3) = e*Sc_x*p / oec^2 * sin(theta) + p / oec * Sc_x_theta;
            S0(2, 3) = e*Sc_y*p / oec^2 * sin(theta) + p / oec * Sc_y_theta;
            S0(3, 3) = e*Sc_z*p / oec^2 * sin(theta) + p / oec * Sc_z_theta;
            S0(4, 3) = -munapi*Sc_Vx_theta;
            S0(5, 3) = -munapi*Sc_Vy_theta;
            S0(6, 3) = munapi*Sc_Vz_theta;
            
            %     dx = 0.001;
            %     [x0_extr2 y0_extr2 z0_extr2 Vx_extr2 Vy_extr2 Vz_extr2] = ...
            %         get_vector_XV( e, p, theta+dx, omega, Omega, i0);
            %     dS = [x0_extr2; y0_extr2; z0_extr2; Vx_extr2; Vy_extr2; Vz_extr2] - ...
            %             [x0_extr; y0_extr; z0_extr; Vx_extr; Vy_extr; Vz_extr];
            %
            %     S0(1, 3) = dS(1) / dx;
            %     S0(2, 3) = dS(2) / dx;
            %     S0(3, 3) = dS(3) / dx;
            %     S0(4, 3) = dS(4) / dx;
            %     S0(5, 3) = dS(5) / dx;
            %     S0(6, 3) = dS(6) / dx;
            
            % Discriminator for omega
            Sc_x_omega = Sc_x_theta;
            Sc_y_omega = Sc_y_theta;
            Sc_z_omega = Sc_z_theta;
            Sc_Vx_omega = Sc_Vx_theta;
            Sc_Vy_omega = Sc_Vy_theta;
            Sc_Vz_omega = Sc_Vz_theta;
            
            S0(1,4) = p / oec * Sc_x_omega;
            S0(2,4) = p / oec * Sc_y_omega;
            S0(3,4) = p / oec * Sc_z_omega;
            S0(4,4) = munapi * (e*sin(omega)*Sc_x_omega - oec*Sc_Vx_omega);
            S0(5,4) = munapi * (e*sin(omega)*Sc_y_omega - oec*Sc_Vy_omega);
            S0(6,4) = munapi * (e*sin(omega)*Sc_z_omega + oec*Sc_Vz_omega);
            
%             S01(i) = S0(1,4);
%             S02(i) = S0(2,4);
%             S03(i) = S0(3,4);
%             S04(i) = S0(4,4);
%             S05(i) = S0(5,4);
%             S06(i) = S0(6,4);
            
            dx = 1e-6;
            [Xforextr2, Vforextr2] = OF.calcXV( e, p, theta, omega+dx, Omega, i0);
            
            dS = [Xforextr2; Vforextr2] - [OF.Xforextr; OF.Vforextr];
            
            S0(1, 4) = dS(1) / dx;
            S0(2, 4) = dS(2) / dx;
            S0(3, 4) = dS(3) / dx;
            S0(4, 4) = dS(4) / dx;
            S0(5, 4) = dS(5) / dx;
            S0(6, 4) = dS(6) / dx;
            
%             S01_ist(i) = S0(1,4);
%             S02_ist(i) = S0(2,4);
%             S03_ist(i) = S0(3,4);
%             S04_ist(i) = S0(4,4);
%             S05_ist(i) = S0(5,4);
%             S06_ist(i) = S0(6,4);
            
            % Discriminator for Omega
            Sc_x_Omega = -cos(u)*sin(Omega) - sin(u)*cos(Omega)*cos(i0);
            Sc_y_Omega = cos(u)*cos(Omega) - sin(u)*sin(Omega)*cos(i0);
%             Sc_z_Omega = 0;
            Sc_Vx_Omega = -sin(u)*sin(Omega) + cos(u)*cos(Omega)*cos(i0);
            Sc_Vy_Omega = sin(u)*cos(Omega) + cos(u)*sin(Omega)*cos(i0);
%             Sc_Vz_Omega = 0;
            
            S0(1,5) = p/oec*Sc_x_Omega;
            S0(2,5) = p/oec*Sc_y_Omega;
            S0(3,5) = 0;
            S0(4,5) = munapi * (e*sin(omega)*Sc_x_Omega - oec*Sc_Vx_Omega);
            S0(5,5) = munapi * (e*sin(omega)*Sc_y_Omega - oec*Sc_Vy_Omega);
            S0(6,5) = 0;
            
            %     dx = 0.001;
            %     [x0_extr2 y0_extr2 z0_extr2 Vx_extr2 Vy_extr2 Vz_extr2] = ...
            %         get_vector_XV( e, p, theta, omega, Omega+dx, i0);
            %     dS = [x0_extr2; y0_extr2; z0_extr2; Vx_extr2; Vy_extr2; Vz_extr2] - ...
            %             [x0_extr; y0_extr; z0_extr; Vx_extr; Vy_extr; Vz_extr];
            %
            %     S0(1, 5) = dS(1) / dx;
            %     S0(2, 5) = dS(2) / dx;
            %     S0(3, 5) = dS(3) / dx;
            %     S0(4, 5) = dS(4) / dx;
            %     S0(5, 5) = dS(5) / dx;
            %     S0(6, 5) = dS(6) / dx;
            
            % Discriminator for i
            Sc_x_i = sin(u)*sin(Omega)*sin(i0);
            Sc_y_i = -sin(u)*cos(Omega)*sin(i0);
            Sc_z_i = sin(u)*cos(i0);
            Sc_Vx_i = -cos(u)*sin(Omega)*sin(i0);
            Sc_Vy_i = cos(u)*cos(Omega)*sin(i0);
            Sc_Vz_i = cos(u)*cos(i0);
            
            S0(1, 6) = p/oec * Sc_x_i;
            S0(2, 6) = p/oec * Sc_y_i;
            S0(3, 6) = p/oec * Sc_z_i;
            S0(4, 6) = munapi * (e*sin(theta)*Sc_x_i - oec*Sc_Vx_i);
            S0(5, 6) = munapi * (e*sin(theta)*Sc_y_i - oec*Sc_Vy_i);
            S0(6, 6) = munapi * (e*sin(theta)*Sc_z_i + oec*Sc_Vz_i);
            
            %     dx = 0.001;
            %     [x0_extr2 y0_extr2 z0_extr2 Vx_extr2 Vy_extr2 Vz_extr2] = ...
            %         get_vector_XV( e, p, theta, omega, Omega, i0+dx);
            %     dS = [x0_extr2; y0_extr2; z0_extr2; Vx_extr2; Vy_extr2; Vz_extr2] - ...
            %             [x0_extr; y0_extr; z0_extr; Vx_extr; Vy_extr; Vz_extr];
            %
            %     S0(1, 6) = dS(1) / dx;
            %     S0(2, 6) = dS(2) / dx;
            %     S0(3, 6) = dS(3) / dx;
            %     S0(4, 6) = dS(4) / dx;
            %     S0(5, 6) = dS(5) / dx;
            %     S0(6, 6) = dS(6) / dx;
            
            OF.S = (OF.c'*S0')';
            OF.Dextr = OF.F*OF.Dest*OF.F' + OF.GDgG;
        end
        
        function Estimate(OF, X, V)
            
            if abs(X) > 0
                dY = [X; V] - [OF.Xforextr; OF.Vforextr];            
            else
                dY = [X; V]*0;
            end
            OF.dY = dY;
            
            t1 = OF.S'/OF.Dn*OF.S;
            t2 = inv(OF.Dextr);
            OF.Dest = inv(t1 + t2);
            
            OF.Xest = OF.Xextr + OF.Dest*OF.S'/OF.Dn*dY;
            OF.calcXVest();
        end
        
        function [X, V] = calcXV(OF, e, p, theta, omega, Omega, i0)
            
            if e > 0
                 theta =  mod_pm_pi(theta);
                 omega =  mod_pm_pi(omega);
            else
                 e = -e;
                 theta =  mod_pm_pi(-theta);
                 omega =  mod_pm_pi(-omega);
            end
            Omega =  mod_pm_pi(Omega);
            i0 = mod_pm_pi(i0);
            
            %Calc vectors x, y, z and Vx, Vy, Vz for orbital elements
            munapi = sqrt(OF.mu_earth / p);
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
        
        function Prediction(OF)
            if norm(OF.Xforest - [OF.Xpred(1); OF.Ypred(1); OF.Zpred(1)]) > 500
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
        
        function FastInit(OF, X, V)
            Xs(1) = OF.Xest(5);
            Xs(2) = OF.Xest(7);
            Xs(3) = OF.Xest(9);
            Xs(4) = OF.Xest(11);
            Xs(5) = OF.Xest(1);
            Xs(6) = OF.Xest(3)*OF.p_mult;
          
            
            options_solve = optimset('Display','off');  % Turn off display for fsolve
            Xs = fsolve(@(Xfs)(fsolve_Kepler(Xfs, X(1), X(2), X(3),...
                V(1), V(2), V(3))), Xs, options_solve);
            
            OF.Xest(1) = Xs(5);
            OF.Xest(3) = Xs(6) / OF.p_mult;
            OF.Xest(5) = Xs(1);
            OF.Xest(7) = Xs(2);
            OF.Xest(9) = Xs(3);
            OF.Xest(11) = Xs(4);
        end
        
        function M = U1(x)
            M = [1      0       0;
                0      cos(x)  sin(x);
                0      -sin(x) cos(x)];
        end
        
        function M = U3(x)
            M = [cos(x)     sin(x)  0;
                -sin(x)    cos(x)  0;
                0          0       1];
        end
        
        function y = mod_pm_pi( x )
            %MOD_PM_PI mod [-pi; pi];
            y = mod(x + pi, 2*pi) - pi;
        end
        
        function Ysolve = fsolve_Kepler(Xs, Xizm, Yizm, Zizm, VXizm, VYizm, VZizm)
            
            %%%%Xsolve = r, u, Omega, i, d_r, d_u
            % Xsolve = theta, omega_p, Omega, i, e, p
            theta = Xs(1);
            omega_p = Xs(2);
            Omega = Xs(3);
            i = Xs(4);
            e = Xs(5);
            p = Xs(6);
            
            [x, y, z, Vx, Vy, Vz] = OF.calcXV( e, p, theta, omega_p, Omega, i);
            
            ErrX = Xizm - x;
            ErrY = Yizm - y;
            ErrZ = Zizm - z;
            
            ErrVx = VXizm - Vx;
            ErrVy = VYizm - Vy;
            ErrVz = VZizm - Vz;
            
            Ysolve = [ErrX, ErrY, ErrZ, ErrVx, ErrVy, ErrVz];
            
        end
        
    end
end

