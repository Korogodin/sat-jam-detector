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
if SR.TimeOfWeek(end) > SR.beginingTime
    xlim(hA, [SR.beginingTime/1000 t(end)]);
end

footer; % DO NOT EDIT
end