clear;
close all
clc;

load Test2.mat
Nmod = length(X);

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
ImitJamAlert = zeros(1, Nmod);

x_izm = nan(1, Nmod);
y_izm = nan(1, Nmod);
z_izm = nan(1, Nmod);
Vx_izm = nan(1, Nmod);
Vy_izm = nan(1, Nmod);
Vz_izm = nan(1, Nmod);

X_ECI = nan(1, Nmod);
Y_ECI = nan(1, Nmod);
Z_ECI = nan(1, Nmod);
Vx_ECI = nan(1, Nmod);
Vy_ECI = nan(1, Nmod);
Vz_ECI = nan(1, Nmod);

theta_j = nan(1, Nmod);
omega_j = nan(1, Nmod);
Omega_j = nan(1, Nmod);
i_j= nan(1, Nmod);
e_j= nan(1, Nmod);
p_j= nan(1, Nmod);

theta_f = nan(1, Nmod);
omega_f = nan(1, Nmod);
Omega_f = nan(1, Nmod);
i_f= nan(1, Nmod);
e_f= nan(1, Nmod);
p_f= nan(1, Nmod);
    
errXestXpred = nan(1, Nmod);
errXizmXpred = nan(1, Nmod);
errXizmXest = nan(1, Nmod);
errXizmXfsolve = nan(1, Nmod);

dXx = 0; dVx = 0;
for j = 1:Nmod
    
    [ECI, V_ECI] = ECEFtoECI(rad2deg((TimeofWeek(j)-TimeofWeek(3000))/1000*w_earth), [X(j); Y(j); Z(j)], [Vx(j); Vy(j); Vz(j)]);
        
    if j > 1
        dTW = TimeofWeek(j) - TimeofWeek(j-1);
    end
    
    X_ECI(j) = ECI(1);
    Y_ECI(j) = ECI(2);
    Z_ECI(j) = ECI(3);
    Vx_ECI(j) = V_ECI(1);
    Vy_ECI(j) = V_ECI(2);
    Vz_ECI(j) = V_ECI(3);
    
    if j > Nmod/10
        dVx = 5;
        dXx = dXx + dVx*T;
    end
    
    x_izm(j) = X_ECI(j) + dXx;
    y_izm(j) = Y_ECI(j);
    z_izm(j) = Z_ECI(j);
    Vx_izm(j) = Vx_ECI(j) + dVx;
    Vy_izm(j) = Vy_ECI(j);
    Vz_izm(j) = Vz_ECI(j);
    
    x_izm_old = x_izm(j);
    y_izm_old = y_izm(j);
    z_izm_old = z_izm(j);

    if (RMS(j) < 10)&&(Solution(j) == 17)
        OF.goodMeasLine = OF.goodMeasLine + 1;
    else
        OF.goodMeasLine = 0;
    end
%     if j == 1
%     if j > 1
%         Osc = [OF.Xest(1); OF.Xest(2)*OF.p_mult; OF.Xest(3); OF.Xest(4); OF.Xest(5); OF.Xest(6)];
%     end
    Osc = ECI2Osc([x_izm(j); y_izm(j); z_izm(j)], [Vx_izm(j); Vy_izm(j); Vz_izm(j)], Osc);
    e_f(j) = Osc(1);
    p_f(j) = Osc(2);
    theta_f(j) = Osc(3);
    omega_f(j) = Osc(4);
    Omega_f(j) = Osc(5);
    i_f(j) = Osc(6);    
%     end
    
    if j == 1
%         OF.FastInit([x0_izm(j); y0_izm(j); z0_izm(j)], [Vx_izm(j); Vy_izm(j); Vz_izm(j)]);
        OF.Xest = [Osc(1); 0; Osc(2) / OF.p_mult; 0; Osc(3); 0; Osc(4); 0; Osc(5); 0; Osc(6); 0];
        OF.calcXVest;
        fprintf('Initialization errors:\n');
        fprintf('\t%.3f m\n', norm([x_izm(j);y_izm(j);z_izm(j)]-OF.Xforest));
        fprintf('\t%.4f m/s\n', norm([Vx_izm(j);Vy_izm(j);Vz_izm(j)]-OF.Vforest));
        fprintf('X = %.3f\n', OF.Xforest(1));
    end
    
    OF.Extrapolate();
    
    if j > 1
        if dTW > 1000
            for j2 = 1:(fix(round(dTW)/1000) - 1)
                OF.Xest = OF.Xextr;
                OF.calcXVest();
                OF.Prediction([x_izm(j); y_izm(j); z_izm(j)]);
                OF.Extrapolate();
            end
        end
    end
    
    if OF.goodMeasLine > 0
        OF.Estimate(Osc);
    else
        OF.Xest = OF.Xextr;
        OF.calcXVest();
        OF.Prediction([x_izm(j); y_izm(j); z_izm(j)]);
    end

    
    errXestXpred(j) =  norm(OF.Xforest - [OF.Xpred(1); OF.Ypred(1); OF.Zpred(1)]);
    errXizmXpred(j) =  norm([x_izm(j); y_izm(j); z_izm(j)] - [OF.Xpred(1); OF.Ypred(1); OF.Zpred(1)]);
    errXizmXest(j) =  norm([x_izm(j); y_izm(j); z_izm(j)] - OF.Xforest);
    errXizmXfsolve(j) =  norm([x_izm(j); y_izm(j); z_izm(j)] - Osc2ECI(Osc));
    
    OF.Prediction([x_izm(j); y_izm(j); z_izm(j)]);
   
    if OF.Xest(1) > 0;
        e_j(j) = OF.Xest(1);
        theta_j(j) = mod_pm_pi(OF.Xest(5));    
        omega_j(j) = mod_pm_pi(OF.Xest(7));
    else
        e_j(j) = -OF.Xest(1);
        theta_j(j) = mod_pm_pi(-OF.Xest(5));    
        omega_j(j) = mod_pm_pi(-OF.Xest(7));
    end    
    p_j(j) = OF.Xest(3)*OF.p_mult;
    Omega_j(j) = mod_pm_pi(OF.Xest(9));
    i_j(j) = mod_pm_pi(OF.Xest(11));

    if j > 500
%         if ImitJamAlert(j-1) == 1
%             ImitJamAlert(j) = 1;
%         else
            ImitJamAlert(j) = OF.ImitJamAlert;
%         end
    else
        ImitJamAlert(j) = 0;
    end
    
    if ~mod(j, fix(Nmod/100))
        fprintf('%.0f %%\n', j/Nmod*100);
    end
   
    figure(1)
    subplot(3,2,1)
    plot(1:j, e_j(1:j), 1:j, e_f(1:j))
    ylabel('e');
    subplot(3,2,2)
    plot(1:j, p_j(1:j), 1:j, p_f(1:j))
    ylabel('p');
    subplot(3,2,3)
    plot(1:j, theta_j(1:j), 1:j, theta_f(1:j))
    ylabel('\theta')
    subplot(3,2,4)
    plot(1:j, omega_j(1:j), 1:j, omega_f(1:j))
    ylabel('\omega')
    subplot(3,2,5)
    plot(1:j, Omega_j(1:j), 1:j, Omega_f(1:j))
    ylabel('\Omega')
    subplot(3,2,6)
    plot(1:j, i_j(1:j), 1:j, i_f(1:j))
    ylabel('i');
    
    figure(2); 
    plot((1:j)*T, ImitJamAlert(1:j));
    xlabel('t, s')
    ylabel('Alert')
    
    figure(3)
    plot((1:j)*T,  errXestXpred(1:j), (1:j)*T,  errXizmXpred(1:j), (1:j)*T,  errXizmXest(1:j), (1:j)*T,  errXizmXfsolve(1:j), (1:j)*T, RMS(1:j))
    xlabel('t, s')
    ylabel('Prediction ans estimation error')
    legend('EstPred', 'IzmPred', 'IzmEst', 'IzmFsolve', 'RMS');
    ylim([0 5000]);
    xlim([0 j*T])

    drawnow
end 





