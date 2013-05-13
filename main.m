try 
    close(handle_fig_main); % Close old output form
end

close all
clear 
clc

globals;

addpath([pwd '/func/tracking']); % Functions for tracking algorithm
addpath([pwd '/func/solve']); % Functions for solving without noise
addpath([pwd '/func/interface']); % Functions for interface

Tmod = 3.2*60*60;  %[s], duration of the simulation

% Magic constants
hF_cont = 0; % Last figure's handles
Font_Size = 8; % Font size for output interface
mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
omega_e = 0.7292115e-4; % [rad/s] Earth's rotation rate
options_solve = optimset('Display','off');  % Turn off display for fsolve

% Load true trajectory of SV
load TrueTrajectory.mat
Nmod = fix(Tmod/T);
if Nmod > Nmod_max
    Nmod = Nmod_max;
    Tmod = (Nmod-1)*T;
end
resize_arrays; % resize arrays for new Tmod

handle_fig_main = fig_main(); % open GUI form
