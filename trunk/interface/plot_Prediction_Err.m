function plot_Prediction_Err(hObject, Kind)
%plot_Prediction_Err Plot discrepancy of prediction and error of estimation

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  [SR.errXizmXpred; SR.errXizmXest; SR.errXizmXfsolve; SR.RMS]);
title(hA, 'Izm-Pred; Izm-Est; Izm-fSolve; RMS');
ylabel(hA, 'norm(Err XYZ), m');
xlabel(hA, 'TimeOfWeek, sec');
ylim(hA, [0 1000]);
grid(hA, 'on');

footer; % DO NOT EDIT
end