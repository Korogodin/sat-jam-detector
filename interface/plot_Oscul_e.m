function plot_Oscul_e(hObject, Kind)
%plot_Oscul_e Plot eccentricity (estimated, fsolved)

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  SR.e_estimated, t,  SR.e_fsolved)
title(hA, 'Eccentricity');
ylabel(hA, 'e');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

footer; % DO NOT EDIT
end