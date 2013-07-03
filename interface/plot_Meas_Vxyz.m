function plot_Meas_Vxyz(hObject, Kind)
%plot_Meas_Vxyz Plot velocities (estimated, measure)

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  [SR.Vx_estimated; SR.Vy_estimated; SR.Vz_estimated], t,  [SR.Vx_measured; SR.Vy_measured; SR.Vz_measured])
title(hA, 'Velocities (ECI)');
ylabel(hA, 'Vx, Vy, Vz, m/s');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

footer; % DO NOT EDIT
end