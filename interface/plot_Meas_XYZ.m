function plot_Meas_XYZ(hObject, Kind)
%plot_Meas_XYZ Plot coordinates (estimated, measure)

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  [SR.x_estimated; SR.y_estimated; SR.z_estimated], t,  [SR.x_measured; SR.y_measured; SR.z_measured])
title(hA, 'Coordinates (ECI)');
ylabel(hA, 'X, Y, Z, m');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

if SR.TimeOfWeek(end) > SR.beginingTime
    xlim(hA, [SR.beginingTime/1000 t(end)]);
end

footer; % DO NOT EDIT
end