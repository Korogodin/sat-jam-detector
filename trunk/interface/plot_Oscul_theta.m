function plot_Oscul_theta(hObject, Kind)
%plot_Oscul_theta Plot true anomaly (estimated, fsolved)

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  SR.theta_estimated, t,  SR.theta_fsolved)
title(hA, 'True anomaly');
ylabel(hA, '\theta, rad');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

footer; % DO NOT EDIT
end