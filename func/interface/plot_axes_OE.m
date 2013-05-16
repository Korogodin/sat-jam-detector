%> ======================================================================
%> @brief plot True trajectory graphs
%> @param handles Main handles struct
%> @param hF Flag (handle) for separate window
%> @param ij Indexies of axes (2-el string)
% ======================================================================
function plot_axes_OE(handles, hF, ij)
globals;

if hF 
    figure(hF); 
    hA = gca;
    XLab = 't, s';
else
    hA = eval(['handles.axes_OE_' ij]);
    set(handles.fig_main,'CurrentAxes', hA)
    XLab = '';
    YLab = '';
end
titl = '';
X1 = tmod; X2 = tmod; X3 = tmod; X4 = tmod;
Y1 = nan(1, Nmod); Y2 = nan(1, Nmod); Y3 = nan(1, Nmod); Y4 = nan(1, Nmod);
leg_str ='';
switch ij
    case '11'
        Y1 = Xist.sqrtA;
        YLab = 'a^{1/2}, m^{1/2}';
        titl = 'Square root of semimajor axis';
    case '12'
        Y1 = Xist.e;
        YLab = 'e';
        titl = 'Eccentricity';
    case '13'
        Y1 = Xist.Crc;
        Y2 = Xist.Crs;
        YLab = 'C_{rc}, C_{rs}, m';        
        titl = 'Harmonic coefficients for radius';
        leg_str = 'legend(hA, ''C_{rc}'', ''C_{rs}'')';
    case '14'
        Y1 = Xist.r;
        Y2 = Xist.A;
        YLab = 'r, a, m';
        titl = 'Full radius, semimajor axis';
        leg_str = 'legend(hA, ''Radius'', ''Semimajor axis'')';
    case '21'
        Y1 = Xist.i0;
        YLab = 'i_0, rad';
        titl = 'Base inclination';
    case '22'
        Y1 = Xist.i_dot;
        YLab = 'i_{dot}, rad/s';
        titl = 'Derivative of the base inclination';
    case '23'
        Y1 = Xist.Cic;
        Y2 = Xist.Cis;
        YLab = 'C_{ic}, C_{is}, rad';        
        titl = 'Harmonic coefficients for inclination';
        leg_str = 'legend(hA, ''C_{ic}'', ''C_{is}'')';        
    case '24'
        Y1 = Xist.i;
        YLab = 'i, rad';
        titl = 'Full inclination';
    case '31'
        Y1 = Xist.M0;
        Y2 = Xist.E;
        Y3 = Xist.theta;
        YLab = 'M_0, E, \theta, rad';        
        titl = 'Mean, eccentricity and true anomaly';
        leg_str = 'legend(hA, ''M_0'', ''E'', ''\theta'')';        
    case '32'
        Y1 = Xist.omega;
        YLab = '\omega, rad';
        titl = 'Argument of periapsis';
    case '33'
        Y1 = Xist.Cuc;
        Y2 = Xist.Cus;
        YLab = 'C_{uc}, C_{us}, rad';         
        titl = 'Harmonic coefficients for argument of latitude';
        leg_str = 'legend(hA, ''C_{uc}'', ''C_{us}'')';
    case '34'
        Y1 = Xist.u;
        YLab = 'u, rad'; 
        titl = 'Argument of latitude';
    case '41'
        Y1 = Xist.Omega;
        YLab = '\Omega, rad';
        titl = 'Longitude of the ascending node';
    case '42'
        Y1 = Xist.Omega_dot;
        YLab = '\Omega_{dot}, rad/s';        
        titl = 'Derivative of the longitude of the ascending node';
    case '43'
        Y1 = Xist.x0;
        Y2 = Xist.y0;
        Y3 = Xist.z0;
        YLab = 'x_0, y_0, z_0, m'; 
        titl = 'Coordinates in the inertial coordinate system';
        leg_str = 'legend(hA, ''x_0'', ''y_0'', ''z_0'')';        
    case '44'
        Y1 = Xist.d_x0;
        Y2 = Xist.d_y0;
        Y3 = Xist.d_z0;
        YLab = 'V_x, V_y, V_z, m'; 
        titl = 'Velocities in the inertial coordinate system';
        leg_str = 'legend(hA, ''V_x'', ''V_y'', ''V_z'')';
end

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
