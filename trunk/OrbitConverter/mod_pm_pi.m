function y = mod_pm_pi( x )
%MOD_PM_PI mod [-pi; pi];
   y = mod(x + pi, 2*pi) - pi;
end