function plot_Prediction_Alert(hObject, Kind)

header; % DO NOT EDIT

global SR

t = SR.TimeOfWeek / 1000;

%             figure(1);
            bar(SR.SNR_GPS)
            ylabel(hA, 'q_{c/no}, dBHz')
            xlabel(hA, 'GPS SV num');
            
% plot(hA, t,  SR.imitJamAlert);
% title(hA, 'Imit Jammer Alert');
% ylabel(hA, 'Alert');
% xlabel(hA, 'TimeOfWeek, sec');
% grid(hA, 'on');
% 
% if SR.TimeOfWeek(end) > SR.beginingTime
%     xlim(hA, [SR.beginingTime/1000 t(end)]);
% end

footer; % DO NOT EDIT
end