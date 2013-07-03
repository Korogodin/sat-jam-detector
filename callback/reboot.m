try
    fclose(instrfind);
end

addpath('../cnavisbinr');

NB = CNavisBinr();
NB.setMode(NB.Mode_Device);
NB.openDevice('/dev/ttyUSB0', 38400);

NB.reset_erase;
NB.reset_erase;