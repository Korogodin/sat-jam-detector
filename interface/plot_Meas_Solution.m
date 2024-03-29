function plot_Meas_Solution(hObject, Kind)
%plot_Meas_Solution Plot solution status

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  SR.Solution)
title(hA, 'Solution status');
ylabel(hA, 'Solution state word');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

if SR.TimeOfWeek(end) > SR.beginingTime
    xlim(hA, [SR.beginingTime/1000 t(end)]);
end

footer; % DO NOT EDIT
end