function NB = controlByMeansDevice( MW )

initSR;

try
    fclose(instrfind);
end

NB = CNavisBinr();
NB.setMode(NB.Mode_Device);

Device = get(SR.hC_Device, 'String');
try
    NB.openDevice(Device, 38400);
catch
    disp('I can not open device');
end
NB.request_88h_bymeansof_27h(1);
pause(0.5);
NB.request_87h_bymeansof_39h(1);
pause(0.5);
NB.sendPacket('B2', [dec2hex(bin2dec('00001000'), 2) dec2hex(bin2dec('00000000'), 2)]); % ECEF
pause(0.5);
NB.request_60h_bymeansof_21h(1);

initControlAlgorithm;

SR.SNR_GPS = NB.SNR_GPS;
SR.SNR_GPS_old = NB.SNR_GPS;
SR.SNR_GPS_oldold = NB.SNR_GPS;

while 1
    ok = NB.getPacketData;
    if ok
        NB.parseData;
        if strcmp(NB.PacketNumber, '88')
            save2SR_NB88;
            eatNewSR;
            MW.replot;
        end
        if strcmp(NB.PacketNumber, '87')
            SR.SNR_GPS_oldold = SR.SNR_GPS_old;
            SR.SNR_GPS_old = SR.SNR_GPS;
            SR.SNR_GPS = NB.SNR_GPS;
            
            if sum((SR.SNR_GPS - SR.SNR_GPS_oldold) < -2) > 3
                if SR.BarrageTime == -1
                    SR.BarrageTime = SR.TimeOfWeek(SR.k88);
                    set(SR.hC_Barrage, 'String', sprintf('Detected at TOW = %.0f', SR.BarrageTime/1000));
                end
            end
        end        
        if strcmp(NB.PacketNumber, '60')
            set(SR.hC_SVNum, 'String', ['SVs: ' num2str(NB.GPS_in_Solution) ' + ' num2str(NB.GLO_in_Solution)]);
        end          
    else
        pause(0.4);
    end
end

end

