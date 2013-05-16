%> ======================================================================
%> @brief Calculate true osculating elements for true trajectory
%> @param handles Main handles struct
%> ======================================================================
function solve_wo_noise( handles )

globals;

% Initial point
Xs(1) = Xist.theta(1);
Xs(2) = Xist.omega(1);
Xs(3) = Xist.Omega(1);
Xs(4) = Xist.i(1);
Xs(5) = Xist.e(1);
Xs(6) = Xist.p(1);
    
set(handles.pb_Solve, 'Enable', 'off');
drawnow
pause(0.01);
for i = 1:Nmod
    
    % Calculating osculating elements by means of the true coordinates and velocities
    Xs = fsolve(@(Xfs)(fsolve_Kepler(Xfs, Xist.x0(i), Xist.y0(i), Xist.z0(i),...
         Xist.d_x0(i), Xist.d_y0(i), Xist.d_z0(i))), Xs, options_solve);
    
    % Saving of results
    Xest_won.e(i) = Xs(5);
    Xest_won.theta(i) = Xs(1);    
    Xest_won.omega(i) = Xs(2);
    Xest_won.Omega(i) = Xs(3);
    Xest_won.i(i) = Xs(4);
    Xest_won.p(i) = Xs(6);
    Xest_won.r(i) = Xest_won.p(i) ./ (1 + Xest_won.e(i)*cos(Xest_won.theta(i)));

    % Get coordinates for oculating elements (for test)
    [Xest_won.x0(i) Xest_won.y0(i) Xest_won.z0(i) ...
        Xest_won.d_x0(i) Xest_won.d_y0(i) Xest_won.d_z0(i)] = ...
        get_vector_XV( Xest_won.e(i), Xest_won.p(i), Xest_won.theta(i), ...
                           Xest_won.omega(i), Xest_won.Omega(i), Xest_won.i(i));    
%     [Xest_won.x0(i) Xest_won.y0(i) Xest_won.z0(i)] = get_vector_XYZ( Xest_won.r(i),...
%                 Xest_won.Omega(i), Xest_won.i(i), Xest_won.theta(i)+Xest_won.omega(i) );

    if ~mod(i, fix(Nmod/100))
        set(handles.txt_Solve, 'String', sprintf('%.0f %%', i/Nmod*100));
    end
    drawnow
    pause(0.01);
end
set(handles.txt_Solve, 'String', sprintf('%.0f %%', 100));

% Calculation of derivatives
Xest_won.d_theta = diff(Xest_won.theta)/dTmod;
Xest_won.d_theta(end+1)=Xest_won.d_theta(end);
Xest_won.d_omega = diff(Xest_won.omega)/dTmod;
Xest_won.d_omega(end+1)=Xest_won.d_omega(end);
Xest_won.d_Omega = diff(Xest_won.Omega)/dTmod;
Xest_won.d_Omega(end+1)=Xest_won.d_Omega(end);
Xest_won.d_i = diff(Xest_won.i)/dTmod;
Xest_won.d_i(end+1)=Xest_won.d_i(end);
Xest_won.d_e = diff(Xest_won.e)/dTmod;
Xest_won.d_e(end+1)=Xest_won.d_e(end);
Xest_won.d_p = diff(Xest_won.p)/dTmod;
Xest_won.d_p(end+1)=Xest_won.d_p(end);
Xest_won.d_r = diff(Xest_won.r)/dTmod;
Xest_won.d_r(end+1)=Xest_won.d_r(end);
    
for i = 1:Nmod
        % Saving of results
    if Xest_won.e(i) > 0;
        Xest_won.e(i) = Xest_won.e(i);
        Xest_won.theta(i) = mod_pm_pi(Xest_won.theta(i));    
        Xest_won.omega(i) = mod_pm_pi(Xest_won.omega(i));
    else
        Xest_won.e(i) = -Xest_won.e(i);
        Xest_won.theta(i) = mod_pm_pi(-Xest_won.theta(i));    
        Xest_won.omega(i) = mod_pm_pi(-Xest_won.omega(i));
    end    
    Xest_won.Omega(i) = mod_pm_pi(Xest_won.Omega(i));
    Xest_won.i(i) = mod_pm_pi(Xest_won.i(i));
end

% Output results to form
for i = 1:5
    for j = 1:5
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'DrawMode', 'fast');
        plot_axes_Track(handles, 0, [num2str(i) num2str(j)]);
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'Tag', ['axes_Track_' num2str(i) num2str(j)]);
        set(eval(['handles.axes_Track_' num2str(i) num2str(j)]), 'ButtonDownFcn', str2func('@(hObject,eventdata)fig_main(''axes_Track_ButtonDownFcn'',hObject,eventdata,guidata(hObject))'));
    end
end

set(handles.pan_Orbital_Elements, 'Visible', 'off');
set(handles.pan_Tracking, 'Visible', 'on');
set(handles.pb_Track, 'Enable', 'on');
set(handles.ed_stdX, 'Enable', 'on');
set(handles.ed_stdV, 'Enable', 'on');

end

