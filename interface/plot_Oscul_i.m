function plot_Oscul_i(hObject, Kind)
%plot_Oscul_i Plot Inclination (estimated, fsolved)

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  SR.i_estimated, t,  SR.i_fsolved)
title(hA, 'Inclination');
ylabel(hA, 'i, rad');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

footer; % DO NOT EDIT
end