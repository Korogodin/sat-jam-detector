function plot_Meas_RMS(hObject, Kind)
%plot_Meas_RMS Plot RMS of measurements 

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  SR.RMS)
title(hA, 'RMS');
ylabel(hA, 'RMS, m');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

footer; % DO NOT EDIT
end