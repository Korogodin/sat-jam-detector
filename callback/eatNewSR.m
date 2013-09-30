if (SR.k88 > SR.beginingk88)&&(SR.beginingk88 ~= -1)

    dTW = SR.TimeOfWeek(SR.k88) - SR.TimeOfWeek(SR.k88-1);
    
    if (SR.RMS(SR.k88) < 10)&&(SR.Solution(SR.k88) == 17)
        OF.goodMeasLine = OF.goodMeasLine + 1;
    else
        OF.goodMeasLine = 0;
    end
        
    X_measured = [SR.x_measured(SR.k88); SR.y_measured(SR.k88); SR.z_measured(SR.k88)];
    V_measured = [SR.Vx_measured(SR.k88); SR.Vy_measured(SR.k88); SR.Vz_measured(SR.k88)];
    Osc = ECI2Osc(X_measured, V_measured, Osc);
    
    if Osc(1) >= 0
        SR.e_fsolved(SR.k88) = Osc(1);
        SR.theta_fsolved(SR.k88) = mod_pm_pi(Osc(3));
        SR.omega_fsolved(SR.k88) = mod_pm_pi(Osc(4));
    else
        SR.e_fsolved(SR.k88) = -Osc(1);
        SR.theta_fsolved(SR.k88) = mod_pm_pi(-Osc(3));
        SR.omega_fsolved(SR.k88) = mod_pm_pi(-Osc(4));
    end
    SR.p_fsolved(SR.k88) = Osc(2);
    SR.Omega_fsolved(SR.k88) = mod_pm_pi(Osc(5));
    SR.i_fsolved(SR.k88) = mod_pm_pi(Osc(6));    

    
    if SR.doFastInit
        SR.doFastInit = 0;
        OF.goodMeasLine = 1;
        OF.Xest = [Osc(1); 0; Osc(2) / OF.p_mult; 0; Osc(3); 0; Osc(4); 0; Osc(5); 0; Osc(6); 0];
        OF.calcXVest;
        fprintf('Initialization errors:\n');
        fprintf('\t%.3f m\n', ...
            norm(X_measured - OF.Xforest));
        fprintf('\t%.4f m/s\n', ...
            norm(V_measured - OF.Vforest));
    end
    
    OF.Extrapolate();
    
    if (dTW > 1000) && abs(dTW < 10000)
        for j2 = 1:(fix(round(dTW)/1000) - 1)
            OF.Xest = OF.Xextr;
            OF.calcXVest();
            OF.Prediction(X_measured);
            OF.Extrapolate();
        end
    end
    
    if OF.goodMeasLine > 0
        OF.Estimate(Osc);
    else
        OF.Xest = OF.Xextr;
        OF.calcXVest();
    end
    
    SR.x_estimated(SR.k88) = OF.Xforest(1);
    SR.y_estimated(SR.k88) = OF.Xforest(2);
    SR.z_estimated(SR.k88) = OF.Xforest(3);

    SR.Vx_estimated(SR.k88) = OF.Vforest(1);
    SR.Vy_estimated(SR.k88) = OF.Vforest(2);
    SR.Vz_estimated(SR.k88) = OF.Vforest(3);

    SR.errXizmXpred(SR.k88) =  norm(X_measured - [OF.Xpred(1); OF.Ypred(1); OF.Zpred(1)]);
    SR.errXizmXest(SR.k88) =  norm(X_measured - OF.Xforest);
    SR.errXizmXfsolve(SR.k88) =  norm(X_measured - Osc2ECI(Osc));
    
    OF.Prediction(X_measured);
   
    if OF.Xest(1) >= 0;
        SR.e_estimated(SR.k88) = OF.Xest(1);
        SR.theta_estimated(SR.k88) = mod_pm_pi(OF.Xest(5));    
        SR.omega_estimated(SR.k88) = mod_pm_pi(OF.Xest(7));
    else
        SR.e_estimated(SR.k88) = -OF.Xest(1);
        SR.theta_estimated(SR.k88) = mod_pm_pi(-OF.Xest(5));    
        SR.omega_estimated(SR.k88) = mod_pm_pi(-OF.Xest(7));
    end    
    SR.p_estimated(SR.k88) = OF.Xest(3)*OF.p_mult;
    SR.Omega_estimated(SR.k88) = mod_pm_pi(OF.Xest(9));
    SR.i_estimated(SR.k88) = mod_pm_pi(OF.Xest(11));

    if (SR.k88 - SR.beginingk88) > 350
            SR.imitJamAlert(SR.k88) = OF.ImitJamAlert;
            if (SR.ImitTime == -1) && (SR.imitJamAlert(SR.k88))
                SR.ImitTime = SR.TimeOfWeek(SR.k88);
                set(SR.hC_Imit, 'String', sprintf('Detected at TOW = %.0f', SR.ImitTime/1000));
            end
    else
        SR.imitJamAlert(SR.k88) = 0;
    end
    
else
    SR.e_fsolved(SR.k88) = NaN;
    SR.theta_fsolved(SR.k88) = NaN;
    SR.omega_fsolved(SR.k88) = NaN;
    SR.p_fsolved(SR.k88) = NaN;
    SR.Omega_fsolved(SR.k88) = NaN;
    SR.i_fsolved(SR.k88) = NaN;    

    SR.e_estimated(SR.k88) = NaN;
    SR.theta_estimated(SR.k88) = NaN;
    SR.omega_estimated(SR.k88) = NaN;
    SR.p_estimated(SR.k88) = NaN;
    SR.Omega_estimated(SR.k88) = NaN;
    SR.i_estimated(SR.k88) = NaN;
    
    SR.x_estimated(SR.k88) = NaN;
    SR.y_estimated(SR.k88) = NaN;
    SR.z_estimated(SR.k88) = NaN;
    SR.Vx_estimated(SR.k88) = NaN;
    SR.Vy_estimated(SR.k88) = NaN;
    SR.Vz_estimated(SR.k88) = NaN;    
    
    SR.errXizmXpred(SR.k88) =  NaN;
    SR.errXizmXest(SR.k88) =  NaN;
    SR.errXizmXfsolve(SR.k88) =  NaN;    
    
    SR.imitJamAlert(SR.k88) = 0;
end
