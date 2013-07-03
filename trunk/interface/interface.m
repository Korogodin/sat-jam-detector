% Controls

Position = [10 10 230 50]; hC = MW.CG.addControl('pushbutton', 'controlByMeansDevice(MW)', 'Use Device', Position);
Position = [10 70 230 30]; SR.hC_Device = MW.CG.addControl('edit', '', '/dev/ttyUSB0', Position);

Position = [10 130 230 50]; hC = MW.CG.addControl('pushbutton', 'controlByMeansDump(MW)', 'Use Dump', Position);
Position = [10 190 230 30]; SR.hC_Dump = MW.CG.addControl('edit', '', 'dump.bin', Position);

Position = [10 250 230 50]; hC = MW.CG.addControl('pushbutton', 'controlByMeansMeasures(MW)', 'Use Measures', Position);
Position = [10 310 230 30]; SR.hC_MeasFile = MW.CG.addControl('edit', '', 'Test.mat', Position);

Position = [10 370 230 30]; SR.hC_Barrage = MW.CG.addControl('text', '', 'Not detected', Position);
Position = [10 410 230 30]; hC = MW.CG.addControl('text', '', 'Barrage Jammer:', Position);

Position = [10 470 230 30]; SR.hC_Imit = MW.CG.addControl('text', '', 'Not detected', Position);
Position = [10 510 230 30]; hC = MW.CG.addControl('text', '', 'Imit Jammer:', Position);

Position = [10 570 230 50]; hC = MW.CG.addControl('pushbutton', 'reinitSignal(MW)', 'Reinit filter', Position);

Position = [10 630 230 30]; SR.hC_RMS = MW.CG.addControl('text', '', 'RMS = ', Position);
Position = [10 670 230 30]; SR.hC_SVNum = MW.CG.addControl('text', '', 'SV Num = ', Position);
% Panels and Axes 


% Panel of osculating parameters
p1 = MW.newPG('Osculating parameters', 'Oscul.param.', 3, 2);
ArrPlace = zeros(3,2); ArrPlace(1, 1) = 1; MW.PG{p1}.newAxes(ArrPlace, 'plot_Oscul_e');
ArrPlace = zeros(3,2); ArrPlace(1, 2) = 1; MW.PG{p1}.newAxes(ArrPlace, 'plot_Oscul_p');
ArrPlace = zeros(3,2); ArrPlace(2, 1) = 1; MW.PG{p1}.newAxes(ArrPlace, 'plot_Oscul_theta');
ArrPlace = zeros(3,2); ArrPlace(2, 2) = 1; MW.PG{p1}.newAxes(ArrPlace, 'plot_Oscul_omega');
ArrPlace = zeros(3,2); ArrPlace(3, 1) = 1; MW.PG{p1}.newAxes(ArrPlace, 'plot_Oscul_Omega');
ArrPlace = zeros(3,2); ArrPlace(3, 2) = 1; MW.PG{p1}.newAxes(ArrPlace, 'plot_Oscul_i');

% Panel of received measurements
p2 = MW.newPG('Measurements', 'Measurements', 2, 2);
ArrPlace = zeros(2,2); ArrPlace(1, 1) = 1; MW.PG{p2}.newAxes(ArrPlace, 'plot_Meas_XYZ');
ArrPlace = zeros(2,2); ArrPlace(1, 2) = 1; MW.PG{p2}.newAxes(ArrPlace, 'plot_Meas_Vxyz');
ArrPlace = zeros(2,2); ArrPlace(2, 1) = 1; MW.PG{p2}.newAxes(ArrPlace, 'plot_Meas_Solution');
ArrPlace = zeros(2,2); ArrPlace(2, 2) = 1; MW.PG{p2}.newAxes(ArrPlace, 'plot_Meas_RMS');

% Panel of predictions
p3 = MW.newPG('Prediction', 'Prediction', 2, 1);
ArrPlace = zeros(2,1); ArrPlace(1, 1) = 1; MW.PG{p3}.newAxes(ArrPlace, 'plot_Prediction_Err');
ArrPlace = zeros(2,1); ArrPlace(2, 1) = 1; MW.PG{p3}.newAxes(ArrPlace, 'plot_Prediction_Alert');
