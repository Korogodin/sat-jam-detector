 
e = nan(1, length(X));
p = nan(1, length(X));
theta = nan(1, length(X));
omega = nan(1, length(X));
Omega = nan(1, length(X));
i = nan(1, length(X));

xc = nan(1, length(X));
yc = nan(1, length(X));
zc = nan(1, length(X));
Vxc = nan(1, length(X));
Vyc = nan(1, length(X));
Vzc = nan(1, length(X));

for j = 1:length(X)
    Osc = ECI2Osc([X(j); Y(j); Z(j)], [Vx(j); Vy(j); Vz(j)], [0; 7e6; 0; 0; 0; 0.85]);
    e(j) = Osc(1);
    p(j) = Osc(2);
    theta(j) = Osc(3);
    omega(j) = Osc(4);
    Omega(j) = Osc(5);
    i(j) = Osc(6);
    
    [Xc, Vc] = Osc2ECI(Osc);
    xc(j) = Xc(1);
    yc(j) = Xc(2);
    zc(j) = Xc(3);
    Vxc(j) = Vc(1);
    Vyc(j) = Vc(2);
    Vzc(j) = Vc(3);
    
    figure(10)
    subplot(3,2,1)
    plot(1:j, e(1:j))
    ylabel('e');
    subplot(3,2,2)
    plot(1:j, p(1:j))
    ylabel('p');
    subplot(3,2,3)
    plot(1:j, theta(1:j))
    ylabel('\theta')
    subplot(3,2,4)
    plot(1:j, omega(1:j))
    ylabel('\omega')
    subplot(3,2,5)
    plot(1:j, Omega(1:j))
    ylabel('\Omega')
    subplot(3,2,6)
    plot(1:j, i(1:j))
    ylabel('i');

    figure(11)
    subplot(3,2,1)
    plot(1:j, xc(1:j) - X(1:j))
    ylabel('\Delta X');
    subplot(3,2,2)
    plot(1:j, yc(1:j) - Y(1:j))
    ylabel('\Delta Y');
    subplot(3,2,3)
    plot(1:j, zc(1:j) - Z(1:j))
    ylabel('\Delta Z')
    subplot(3,2,4)
    plot(1:j, Vxc(1:j) - Vx(1:j))
    ylabel('\Delta Vx')
    subplot(3,2,5)
    plot(1:j, Vyc(1:j) - Vy(1:j))
    ylabel('\Delta Vy')
    subplot(3,2,6)
    plot(1:j, Vzc(1:j) - Vz(1:j))
    ylabel('\Delta Vz');
    drawnow
end
