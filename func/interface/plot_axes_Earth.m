%> ======================================================================
%> @brief 3D-Earth, SV's orbits
%> @param handles Main handles struct
%> @param hF Flag (handle) for separate window
% ======================================================================
function plot_axes_Earth(handles, hF)
globals;

R_earth_pol = 6356777; % Polar radius of Earth
R_earth_equa = 6378160; % Equatorial radius of Earth
Earth_axe_angl = deg2rad(23); 

[x,y,z] = sphere(50);
x = R_earth_equa * x; y = R_earth_equa * y; z = R_earth_pol * z;
load('topo.mat','topo','topomap1');
% cla reset
% axis square off
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
alpha_earth = pi/1;
x_new = x*cos(alpha_earth) + y*sin(alpha_earth);
y_new = -x*sin(alpha_earth) + y*cos(alpha_earth);

if hF 
    figure(hF); 
    hA = gca;
else
    set(handles.fig_main,'CurrentAxes', handles.axes_3D)
    hA = handles.axes_3D;
end

cla(hA);
surface(x_new, y_new, z, props);
light('position',[-1 -1 1*tan(Earth_axe_angl)]*R_earth_equa*2);
% campos([+1 0 -1*tan(Earth_axe_angl)]*R_earth_equa*2);
% camtarget([0 0 0])
view(3)        
view(40,-8);
hold(hA, 'on');
plot3(hA, Xest.x, Xest.y, Xest.z, Xist.x, Xist.y, Xist.z, 'LineWidth', 2);   
hold(hA, 'off');

xlim([-8e6 8e6]); ylim([-8e6 8e6]); zlim([-8e6 8e6]);
title(hA, 'Space View');
if ~hF
    set(hA, 'XTick', []);
    set(hA, 'YTick', []);
    set(hA, 'ZTick', []);
    set(hA, 'FontSize', Font_Size);
else
    xlabel(hA, 'x, m');
    ylabel(hA, 'y, m');
    zlabel(hA, 'z, m');
    grid(hA, 'on');
end