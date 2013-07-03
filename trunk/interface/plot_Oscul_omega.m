function plot_Oscul_omega(hObject, Kind)
%plot_Oscul_omega Plot Argument of periapsis (estimated, fsolved)

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  SR.omega_estimated, t,  SR.omega_fsolved)
title(hA, 'Argument of periapsis');
ylabel(hA, '\omega, rad');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

footer; % DO NOT EDIT
end