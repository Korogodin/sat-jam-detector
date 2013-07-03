    [ECI, V_ECI] = ECEFtoECI(rad2deg((TimeofWeek(j)-TimeofWeek(1))/1000*w_earth), [X(j); Y(j); Z(j)], [Vx(j); Vy(j); Vz(j)]);
        
    if j > 1
        dTW = TimeofWeek(j) - TimeofWeek(j-1);
    end
    
    X_ECI(j) = ECI(1);
    Y_ECI(j) = ECI(2);
    Z_ECI(j) = ECI(3);
    Vx_ECI(j) = V_ECI(1);
    Vy_ECI(j) = V_ECI(2);
    Vz_ECI(j) = V_ECI(3);
    
    
    x_izm(j) = X_ECI(j) ;
    y_izm(j) = Y_ECI(j);
    z_izm(j) = Z_ECI(j);
    Vx_izm(j) = Vx_ECI(j) ;
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
    
  
%     figure(1)
%     subplot(3,2,1)
%     plot(1:j, e_j(1:j), 1:j, e_f(1:j))
%     ylabel('e');
%     subplot(3,2,2)
%     plot(1:j, p_j(1:j), 1:j, p_f(1:j))
%     ylabel('p');
%     subplot(3,2,3)
%     plot(1:j, theta_j(1:j), 1:j, theta_f(1:j))
%     ylabel('\theta')
%     subplot(3,2,4)
%     plot(1:j, omega_j(1:j), 1:j, omega_f(1:j))
%     ylabel('\omega')
%     subplot(3,2,5)
%     plot(1:j, Omega_j(1:j), 1:j, Omega_f(1:j))
%     ylabel('\Omega')
%     subplot(3,2,6)
%     plot(1:j, i_j(1:j), 1:j, i_f(1:j))
%     ylabel('i');
%     
%     figure(2); 
%     plot((1:j)*T, ImitJamAlert(1:j));
%     xlabel('t, s')
%     ylabel('Alert')
%     
%     figure(3)
%     plot((1:j)*T,  errXestXpred(1:j), (1:j)*T,  errXizmXpred(1:j), (1:j)*T,  errXizmXest(1:j), (1:j)*T,  errXizmXfsolve(1:j), (1:j)*T, RMS(1:j))
%     xlabel('t, s')
%     ylabel('Prediction ans estimation error')
%     legend('EstPred', 'IzmPred', 'IzmEst', 'IzmFsolve', 'RMS');
%     ylim([0 5000]);
%     xlim([0 j*T])
