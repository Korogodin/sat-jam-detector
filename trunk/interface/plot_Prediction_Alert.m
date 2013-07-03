function plot_Prediction_Alert(hObject, Kind)

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  SR.imitJamAlert);
title(hA, 'Imit Jammer Alert');
ylabel(hA, 'Alert');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

footer; % DO NOT EDIT
end