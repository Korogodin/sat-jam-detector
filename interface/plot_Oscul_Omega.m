function plot_Oscul_Omega(hObject, Kind)
%plot_Oscul_Omega Plot Longitude of the ascending node (estimated, fsolved)

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  SR.Omega_estimated, t,  SR.Omega_fsolved)
title(hA, 'Longitude of the ascending node');
ylabel(hA, '\Omega, rad');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

if SR.TimeOfWeek(end) > SR.beginingTime
    xlim(hA, [SR.beginingTime/1000 t(end)]);
end

footer; % DO NOT EDIT
end