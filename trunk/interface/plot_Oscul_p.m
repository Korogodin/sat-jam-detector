function plot_Oscul_p(hObject, Kind)
%plot_Oscul_p Plot focal parameter (estimated, fsolved)

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;
plot(hA, t,  SR.p_estimated, t,  SR.p_fsolved)
title(hA, 'Focal parameter');
ylabel(hA, 'p, m');
xlabel(hA, 'TimeOfWeek, sec');
grid(hA, 'on');

if SR.TimeOfWeek(end) > SR.beginingTime
    xlim(hA, [SR.beginingTime/1000 t(end)]);
end

footer; % DO NOT EDIT
end