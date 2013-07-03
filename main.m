try 
    close(MW.handle); % Close old output form
end

global SR
try
    clear SR
end

close all
clear
clc

rng('shuffle'); % Reinit for randn

addpath([pwd '/OrbitConverter']); 
addpath([pwd '/callback']); % Callbacks for Control Panel
addpath([pwd '/basic-interface']); % Functions for base interface
addpath([pwd '/interface']); % Functions for interface
addpath([pwd '/cnavisbinr']); % Class of NAVIS binary protocol 

MW = CMainWindow('SV GNSS Control');
if MW.handle == 0
    clear MW;
    return;
end

global SR
interface;