
fid = fopen('LO.umt', 'w');
% fprintf(fid, 'Time_ms,Pos_X,Pos_Y,Pos_Z,Vel_X,Vel_Y,Vel_Z\n');

Xist.d_x = diff(Xist.x); Xist.d_x(end+1) = Xist.d_x(end);
Xist.d_y = diff(Xist.y); Xist.d_y(end+1) = Xist.d_y(end);
Xist.d_z = diff(Xist.z); Xist.d_z(end+1) = Xist.d_z(end);

Xist.dd_x = diff(Xist.d_x); Xist.dd_x(end+1) = Xist.dd_x(end);
Xist.dd_y = diff(Xist.d_y); Xist.dd_y(end+1) = Xist.dd_y(end);
Xist.dd_z = diff(Xist.d_z); Xist.dd_z(end+1) = Xist.dd_z(end);

for i = 1:length(tmod)
    fprintf(fid, datestr(i/24/60/60,  'HH:MM:SS'));
    fprintf(fid, ',MOT,v1_m1');
%     fprintf(fid, '%.0f', tmod(i)*1000);
    fprintf(fid, ',%.4f', Xist.x(i));
    fprintf(fid, ',%.4f', Xist.y(i));
    fprintf(fid, ',%.4f', Xist.z(i));
    fprintf(fid, ',%.4f', Xist.d_x(i));
    fprintf(fid, ',%.4f', Xist.d_y(i));
    fprintf(fid, ',%.4f', Xist.d_z(i));    
    fprintf(fid, ',%.4f', Xist.dd_x(i));
    fprintf(fid, ',%.4f', Xist.dd_y(i));
    fprintf(fid, ',%.4f', Xist.dd_z(i));    
%     fprintf(fid, ',0,0,0,0,0,0,0,0,0');
    fprintf(fid, '\n');
end