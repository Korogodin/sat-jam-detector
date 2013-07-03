function NB = controlByMeansDump( MW )

initSR;

NB = CNavisBinr();
NB.setMode(NB.Mode_Dump);

Dump = get(SR.hC_Dump, 'String');
NB.openDump(Dump);

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
    end
end

end

