global SR

SR.k88 = SR.k88 + 1;

[ECI, V_ECI] = ECEFtoECI(rad2deg(NB.TimeOfWeek/1000*w_earth), [NB.X; NB.Y; NB.Z], [NB.Vx; NB.Vy; NB.Vz]);

SR.x_measured(SR.k88) = ECI(1);
SR.y_measured(SR.k88) = ECI(2);
SR.z_measured(SR.k88) = ECI(3);

SR.Vx_measured(SR.k88) = V_ECI(1);
SR.Vy_measured(SR.k88) = V_ECI(2);
SR.Vz_measured(SR.k88) = V_ECI(3);

SR.RMS(SR.k88) = NB.RMS;
set(SR.hC_RMS, 'String', ['RMS = ' num2str(NB.RMS) ' m']);
SR.Solution(SR.k88) = NB.Solution;

SR.TimeOfWeek(SR.k88) = NB.TimeOfWeek;

if (NB.Solution==17)&&(SR.beginingTime == -1)
    SR.beginingTime = NB.TimeOfWeek;
    SR.beginingk88 = SR.k88;
    SR.doFastInit = 1;
end


