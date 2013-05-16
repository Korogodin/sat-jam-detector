%> ======================================================================
%> @brief Tracking alghorithm
%> @param handles Main handles struct
%> ======================================================================
function track_with_noise( handles, std_x, std_V )

set(handles.pb_Track, 'Enable', 'off'); 
pause(0.01);

globals;

% Filter step
T = dTmod;

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

Xextr.X =  [Xest_won.e(2); 0; Xest_won.p(2)/p_mult; 0; Xest_won.theta(2); 0; Xest_won.omega(2); 0; Xest_won.Omega(2); 0; Xest_won.i(2); 0];
% Xextr.X =  [0.01; 0; Xest_won.p(1)/p_mult; 0; 1; Xest_won.d_theta(1)*1.1; 0; 0; 0; Xest_won.d_Omega(1)*0.9; Xest_won.i(1)*0.9; Xest_won.d_i(1)];
% Xextr.X =  [0.01; 0; 7e6/p_mult; 0; 1; Xest_won.d_theta(1)*1.1; 0; 0; 0; Xest_won.d_Omega(1)*0.9; Xest_won.i(1)*0.9; Xest_won.d_i(1)];
Xest.X = Xextr.X * 1.5;


% RMS of shaping noises
std_e = 5e-7 / 15*dTmod;
std_p = 5 / 15*dTmod / p_mult; % [m]
std_theta = 1e-5 / 15*dTmod; % [rad]
std_omega = 1e-6 / 15*dTmod; % [rad]
std_Omega = 1e-8 / 15*dTmod; % [rad]
std_i = 1e-8 / 15*dTmod; % [rad]

Dest = [std_e^2*1e1     0           0           0               0               0               0               0               0               0               0           0
            0           std_e^2*1e2 0           0               0               0               0               0               0               0               0           0
            0           0           std_p^2*1e3 0               0               0               0               0               0               0               0           0
            0           0           0           std_p^2*1e4     0               0               0               0               0               0               0           0
            0           0           0           0               std_theta^2*1e2 0               0               0               0               0               0           0
            0           0           0           0               0               std_theta^2*1e0 0               0               0               0               0           0
            0           0           0           0               0               0               std_omega^2*1e2 0               0               0               0           0
            0           0           0           0               0               0               0               std_omega^2*1e2 0               0               0           0
            0           0           0           0               0               0               0               0               std_Omega^2*1e2 0               0           0
            0           0           0           0               0               0               0               0               0               std_Omega^2*1e2 0           0
            0           0           0           0               0               0               0               0               0               0               std_i^2*1e2 0
            0           0           0           0               0               0               0               0               0               0               0           std_i^2*1e2];

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

% std_x = 10 / sqrt(dTmod)  ;
% std_V = 0.01 / sqrt(dTmod);

Dn = [std_x^2  0         0          0       0       0
      0        std_x^2   0          0       0       0
      0        0         std_x^2    0       0       0
      0        0         0          std_V^2 0       0
      0        0         0          0       std_V^2 0
      0        0         0          0       0       std_V^2];
 
c = [1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0 0; 
     0 0 0 0 0 0 0 0 0 0 1 0];
 
OF = COsculFilter(T, Xest.X, Dest, p_mult, Dg, Dn, F, G, c);
global ImitJamAlert
ImitJamAlert = zeros(1, Nmod);
dXx = 0;
for i = 1:Nmod
    
    x0_izm = Xist.x0(i) + std_x * randn(1,1);
    y0_izm = Xist.y0(i) + std_x * randn(1,1);
    z0_izm = Xist.z0(i) + std_x * randn(1,1);
    Vx_izm = Xist.d_x0(i) + std_V * randn(1,1);    
    Vy_izm = Xist.d_y0(i) + std_V * randn(1,1);
    Vz_izm = Xist.d_z0(i) + std_V * randn(1,1);

    if i == 1
        OF.FastInit([x0_izm; y0_izm; z0_izm], [Vx_izm; Vy_izm; Vz_izm]);
    end
    
    OF.Extrapolate();
    if i > Nmod/2
        dVx = 5;
        dXx = dXx + dVx*T;
        OF.Estimate([x0_izm+dXx; y0_izm; z0_izm], [Vx_izm+dVx; Vy_izm; Vz_izm]);
    else
        OF.Estimate([x0_izm; y0_izm; z0_izm], [Vx_izm; Vy_izm; Vz_izm]);
    end
    OF.Prediction();
    
    if Xest.X(1) > 0;
        Xest.e(i) = OF.Xest(1);
        Xest.theta(i) = mod_pm_pi(OF.Xest(5));    
        Xest.omega(i) = mod_pm_pi(OF.Xest(7));
    else
        Xest.e(i) = -OF.Xest(1);
        Xest.theta(i) = mod_pm_pi(-OF.Xest(5));    
        Xest.omega(i) = mod_pm_pi(-OF.Xest(7));
    end    
    Xest.p(i) = OF.Xest(3)*OF.p_mult;
    Xest.Omega(i) = mod_pm_pi(OF.Xest(9));
    Xest.i(i) = mod_pm_pi(OF.Xest(11));
    
    [Xest.x0(i) Xest.y0(i) Xest.z0(i) Xest.d_x0(i) Xest.d_y0(i) Xest.d_z0(i)] = ...
        get_vector_XV( Xest.e(i), Xest.p(i), Xest.theta(i), Xest.omega(i), Xest.Omega(i), Xest.i(i) );
    [Xest.x(i) Xest.y(i) Xest.z(i) tmp1 tmp2 tmp3] = ...
        get_vector_XV( Xest.e(i), Xest.p(i), Xest.theta(i), Xest.omega(i), Xest.Omega(i) - omega_e*tmod(i), Xest.i(i) );
    
    Xest.d_e(i) = OF.Xest(2);
    Xest.d_p(i) = OF.Xest(4)*OF.p_mult;
    Xest.d_theta(i) = OF.Xest(6);
    Xest.d_omega(i) = OF.Xest(8);
    Xest.d_Omega(i) = OF.Xest(10);
    Xest.d_i(i) = OF.Xest(12);
    
    if i > 70
        if ImitJamAlert(i-1) == 1
            ImitJamAlert(i) = 1;
        else
            ImitJamAlert(i) = OF.ImitJamAlert;
        end
    else
        ImitJamAlert(i) = 0;
    end
    
    if ~mod(i, fix(Nmod/100))
        set(handles.txt_Track, 'String', sprintf('%.0f %%', i/Nmod*100));
    end
    drawnow
    pause(0.01);    
end 
set(handles.txt_Track, 'String', sprintf('%.0f %%', 100));

set(handles.pb_Track, 'Enable', 'on'); 

set(handles.pan_Orbital_Elements, 'Visible', 'off');
set(handles.pan_Tracking, 'Visible', 'on');

% Output results to form
for i = 1:5
    for j = 1:5
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'DrawMode', 'fast');
        plot_axes_Track(handles, 0, [num2str(i) num2str(j)]);
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'Tag', ['axes_Track_' num2str(i) num2str(j)]);
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'ButtonDownFcn', str2func('@(hObject,eventdata)fig_main(''axes_Track_ButtonDownFcn'',hObject,eventdata,guidata(hObject))'));
    end
end
plot_axes_Earth(handles, 0);
draw_Errors(handles);


figure(1); 
plot((1:Nmod)*T, ImitJamAlert);
xlabel('t, s')
ylabel('Alert')
end


