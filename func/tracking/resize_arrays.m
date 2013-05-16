
tmod = tmod(1:Nmod);

% True vector
Xist.dn = Xist.dn(1:Nmod);
Xist.M0 = Xist.M0(1:Nmod);
Xist.e = Xist.e(1:Nmod);
Xist.d_e = Xist.d_e(1:Nmod);
Xist.dd_e = Xist.dd_e(1:Nmod);
Xist.A = Xist.A(1:Nmod);
Xist.d_A = Xist.d_A(1:Nmod);
Xist.dd_A = Xist.dd_A(1:Nmod);
Xist.sqrtA = Xist.sqrtA(1:Nmod);
Xist.Omega = Xist.Omega(1:Nmod);
Xist.d_Omega = Xist.d_Omega(1:Nmod);
Xist.dd_Omega = Xist.dd_Omega(1:Nmod);
Xist.lambda = Xist.lambda(1:Nmod);           
Xist.d_lambda = Xist.d_lambda(1:Nmod);
Xist.dd_lambda = Xist.dd_lambda(1:Nmod);
Xist.Omega_dot = Xist.Omega_dot(1:Nmod);
Xist.Cis = Xist.Cis(1:Nmod);
Xist.Cic = Xist.Cic(1:Nmod);
Xist.Crs = Xist.Crs(1:Nmod);
Xist.Crc = Xist.Crc(1:Nmod);
Xist.Cus = Xist.Cus(1:Nmod);
Xist.Cuc = Xist.Cuc(1:Nmod);
Xist.omega = Xist.omega(1:Nmod);
Xist.d_omega = Xist.d_omega(1:Nmod);
Xist.dd_omega = Xist.dd_omega(1:Nmod);
Xist.i0 = Xist.i0(1:Nmod);           
Xist.i_dot = Xist.i_dot(1:Nmod);
Xist.i = Xist.i(1:Nmod);
Xist.d_i = Xist.d_i(1:Nmod);
Xist.dd_i = Xist.dd_i(1:Nmod);
Xist.u = Xist.u(1:Nmod);
Xist.d_u = Xist.d_u(1:Nmod);
Xist.dd_u = Xist.dd_u(1:Nmod);
Xist.r = Xist.r(1:Nmod);
Xist.d_r = Xist.d_r(1:Nmod);
Xist.dd_r = Xist.dd_r(1:Nmod);
Xist.E = Xist.E(1:Nmod);
Xist.x0 = Xist.x0(1:Nmod);
Xist.y0 = Xist.y0(1:Nmod);
Xist.z0 = Xist.z0(1:Nmod);
Xist.theta1 = Xist.theta1(1:Nmod);
Xist.theta = Xist.theta(1:Nmod);
Xist.d_theta = Xist.d_theta(1:Nmod);
Xist.dd_theta = Xist.dd_theta(1:Nmod);
Xist.p = Xist.p(1:Nmod);
Xist.d_p = Xist.d_p(1:Nmod);
Xist.dd_p = Xist.dd_p(1:Nmod);
Xist.d_x0 = Xist.d_x0(1:Nmod);
Xist.d_y0 = Xist.d_y0(1:Nmod);
Xist.d_z0 = Xist.d_z0(1:Nmod);
Xist.x = Xist.x(1:Nmod);
Xist.y = Xist.y(1:Nmod);
Xist.z = Xist.z(1:Nmod);


% Estimations and extrapolations of Kalman 
Xest.e = nan(1, Nmod);
Xest.d_e = nan(1, Nmod);
Xest.p = nan(1, Nmod);
Xest.d_p = nan(1, Nmod);
Xest.theta = nan(1, Nmod);
Xest.d_theta = nan(1, Nmod);
Xest.omega = nan(1, Nmod);
Xest.d_omega = nan(1, Nmod);
Xest.Omega = nan(1, Nmod);
Xest.d_Omega = nan(1, Nmod);
Xest.i = nan(1, Nmod);
Xest.d_i = nan(1, Nmod);
Xest.x0 = nan(1, Nmod);
Xest.y0 = nan(1, Nmod);
Xest.z0 = nan(1, Nmod);
Xest.d_x0 = nan(1, Nmod);
Xest.d_y0 = nan(1, Nmod);
Xest.d_z0 = nan(1, Nmod);
Xest.lambda = nan(1, Nmod);
Xest.x = nan(1, Nmod);
Xest.y = nan(1, Nmod);
Xest.z = nan(1, Nmod);
Xest.X = nan(12,1);
Xextr = Xest; % - extrapolation

% Without noise solution
Xest_won = Xest; 
