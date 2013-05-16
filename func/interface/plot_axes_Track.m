%> ======================================================================
%> @brief plot Tracking graphs
%> @param handles Main handles struct
%> @param hF Flag (handle) for separate window
%> @param ij Indexies of axes (2-el string)
% ======================================================================
function plot_axes_Track(handles, hF, ij)
globals;

if hF 
    figure(hF); 
    hA = gca;
    XLab = 't, s';
else
    hA = eval(['handles.axes_Track_' ij]);
    set(handles.fig_main,'CurrentAxes', hA)
    XLab = '';
    YLab = '';
end

titl = '';
X1 = tmod; X2 = tmod; X3 = tmod; X4 = tmod; X5 = tmod; X6 = tmod;
Y1 = nan(1, Nmod); Y2 = nan(1, Nmod); Y3 = nan(1, Nmod); Y4 = nan(1, Nmod); Y5 = nan(1, Nmod); Y6 = nan(1, Nmod);
switch ij
    case '11'
        Y1 = Xest.e;
        Y2 = Xest_won.e;
        YLab = 'e';
        titl = 'Eccentricity';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '12'
        Y1 = Xest.d_e;
        Y2 = Xest_won.d_e;
        YLab = 'e''';
        titl = 'Eccentricity rate';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '13'
        Y1 = Xest.p;
        Y2 = Xest_won.p;
        YLab = 'p, m';
        titl = 'Focal parameter';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '14'
        Y1 = Xest.d_p;
        Y2 = Xest_won.d_p;
        YLab = 'p'', m/s';
        titl = 'Focal parameter rate';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '21'
        Y1 = Xest.theta;
        Y2 = Xest_won.theta;
        YLab = '\theta, rad';
        titl = 'True anomaly';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '22'
        Y1 = Xest.d_theta;
        Y2 = Xest_won.d_theta;
        YLab = '\theta'', rad/s';
        titl = 'True anomaly rate';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '23'
        Y1 = Xest.omega;
        Y2 = Xest_won.omega;
        YLab = '\omega, rad';
        titl = 'Argument of periapsis';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '24'
        Y1 = Xest.d_omega;
        Y2 = Xest_won.d_omega;
        YLab = '\omega'', rad/s';
        titl = 'Argument of periapsis rate';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '31'
        Y1 = Xest.Omega;
        Y2 = Xest_won.Omega;
        YLab = '\Omega, rad';
        titl = 'Longitude of the ascending node';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '32'
        Y1 = Xest.d_Omega;
        Y2 = Xest_won.d_Omega;
        YLab = '\Omega'', rad/s';
        titl = 'Longitude of the ascending node rate';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '33'
        Y1 = Xest.i;
        Y2 = Xest_won.i;
        YLab = 'i, rad';
        titl = 'Inclination';
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '34'
        Y1 = Xest.d_i;
        Y2 = Xest_won.d_i;
        YLab = 'i'', rad/s';
        titl = 'Inclination rate'; 
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '41'
        Y1 = Xest.x0;
        Y2 = Xest.y0;
        Y3 = Xest.z0;
        Y4 = Xest_won.x0;
        Y5 = Xest_won.y0;
        Y6 = Xest_won.z0;
        YLab = 'x_0, y_0, z_0, m';
        titl = 'Coordinates in the inertial coordinate system';        
        leg_str = 'legend(hA, ''Estimation x_0'', ''Estimation y_0'', ''Estimation z_0'', ''True x_0'', ''True y_0'', ''True z_0'')';
    case '42'
        Y1 = Xest.d_x0;
        Y2 = Xest.d_y0;
        Y3 = Xest.d_z0;
        Y4 = Xest_won.d_x0;
        Y5 = Xest_won.d_y0;
        Y6 = Xest_won.d_z0;
        YLab = 'V_x, V_y, V_z, m/s';
        titl = 'Velocities in the inertial coordinate system';        
        leg_str = 'legend(hA, ''Estimation V_x'', ''Estimation V_y'', ''Estimation V_z'', ''True V_x'', ''True V_y'', ''True V_z'')';
    case '43'
        Y1 = sqrt((Xest.x0 - Xist.x0).^2 + (Xest.y0 - Xist.y0).^2 + (Xest.z0 - Xist.z0).^2);
        Y2 = sqrt((Xest_won.x0 - Xist.x0).^2 + (Xest_won.y0 - Xist.y0).^2 + (Xest_won.z0 - Xist.z0).^2);
        YLab = 'Error XYZ, m';
        titl = 'Error of coordinate estimation';     
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
    case '44'
        Y1 = sqrt((Xest.d_x0 - Xist.d_x0).^2 + (Xest.d_y0 - Xist.d_y0).^2 + (Xest.d_z0 - Xist.d_z0).^2);
        Y2 = sqrt((Xest_won.d_x0 - Xist.d_x0).^2 + (Xest_won.d_y0 - Xist.d_y0).^2 + (Xest_won.d_z0 - Xist.d_z0).^2);
        YLab = 'Error V, m/s';
        titl = 'Error of velocity estimation';       
        leg_str = 'legend(hA, ''Estimation'', ''True'')';
end

if isnan(Y6(1))
    if isnan(Y5(1))
        if isnan(Y4(1))
            if isnan(Y3(1))
                if isnan(Y2(1))
                    if isnan(Y1(1))
                        set(hA, 'Visible', 'off');
                    else
                        plot(hA, X1, Y1);
                        set(hA, 'Visible', 'on');
                    end
                else
                    plot(hA, X1, Y1, X2, Y2);
                    set(hA, 'Visible', 'on');
                end
            else
                plot(hA, X1, Y1, X2, Y2, X3, Y3);
                set(hA, 'Visible', 'on');
            end
        else
            plot(hA, X1, Y1, X2, Y2, X3, Y3, X4, Y4);
            set(hA, 'Visible', 'on');
        end
    else
        plot(hA, X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5);
        set(hA, 'Visible', 'on');
    end
else
    plot(hA, X1, Y1, X2, Y2, X3, Y3, X4, Y4, X5, Y5, X6, Y6);
    set(hA, 'Visible', 'on');
end


if hF
    title(hA, titl);
    eval(leg_str);
end

ly = ylabel(YLab);
lx = xlabel(XLab);
grid(hA, 'on');

if ~hF
    set(hA, 'XTick', []);
    set(hA, 'YTick', []);
    set(ly, 'FontSize', Font_Size);
end
