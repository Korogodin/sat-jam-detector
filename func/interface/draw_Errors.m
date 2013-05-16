%> ======================================================================
%> @brief plot Caclulate and draw RMSes
%> @param handles Main handles struct
% ======================================================================
function draw_Errors(handles)
globals;
begin_i = Nmod - 300; end_i = Nmod;

set(handles.txt_Err_11, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.e(begin_i:end_i) - Xist.e(begin_i:end_i)).^2)))));
set(handles.txt_Err_12, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.d_e(begin_i:end_i) - Xist.d_e(begin_i:end_i)).^2)))));
set(handles.txt_Err_21, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.p(begin_i:end_i) - Xist.p(begin_i:end_i)).^2)))));
set(handles.txt_Err_22, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.d_p(begin_i:end_i) - Xist.d_p(begin_i:end_i)).^2)))));
set(handles.txt_Err_31, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.theta(begin_i:end_i) - Xist.theta(begin_i:end_i)).^2)))));
set(handles.txt_Err_32, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.d_theta(begin_i:end_i) - Xist.d_theta(begin_i:end_i)).^2)))));
set(handles.txt_Err_41, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.omega(begin_i:end_i) - Xist.omega(begin_i:end_i)).^2)))));
set(handles.txt_Err_42, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.d_omega(begin_i:end_i) - Xist.d_omega(begin_i:end_i)).^2)))));
set(handles.txt_Err_51, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.Omega(begin_i:end_i) - Xist.Omega(begin_i:end_i)).^2)))));
set(handles.txt_Err_52, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.d_Omega(begin_i:end_i) - Xist.d_Omega(begin_i:end_i)).^2)))));
set(handles.txt_Err_61, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.i(begin_i:end_i) - Xist.i(begin_i:end_i)).^2)))));
set(handles.txt_Err_62, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.d_i(begin_i:end_i) - Xist.d_i(begin_i:end_i)).^2)))));
set(handles.txt_Err_71, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.x0(begin_i:end_i) - Xist.x0(begin_i:end_i)).^2)))));
set(handles.txt_Err_72, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.d_x0(begin_i:end_i) - Xist.d_x0(begin_i:end_i)).^2)))));
set(handles.txt_Err_81, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.y0(begin_i:end_i) - Xist.y0(begin_i:end_i)).^2)))));
set(handles.txt_Err_82, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.d_y0(begin_i:end_i) - Xist.d_y0(begin_i:end_i)).^2)))));
set(handles.txt_Err_91, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.z0(begin_i:end_i) - Xist.z0(begin_i:end_i)).^2)))));
set(handles.txt_Err_92, 'String', num2str(sprintf('%.1E', mean(sqrt((Xest.d_z0(begin_i:end_i) - Xist.d_z0(begin_i:end_i)).^2)))));

end

