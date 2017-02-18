clc
clear
options = optimset('Display','off');

F = 1e3;
R0 = 8.314;
Tc = 1500;
M = 14.4e-3;
k = 1.4;
pa = 1e5;
pe = 3e5;
pc = 100e5;

A2diam = @(x) 2 * sqrt(x/pi);

epsilon = sqrt((k - 1)/2) * (2 / (k+1))^((k+1) / 2 / (k-1)) / (pe/pc)^(1/k) / sqrt(1 - (pe/pc)^((k - 1)/k))
Ve = sqrt( 2*k*R0*Tc / (k-1) / M * (1 - (pe/pc)^((k-1)/k)))
K = sqrt(k) * (2 / (k+1))^((k+1) / 2 / (k-1));
Cstar = sqrt(R0 / M * Tc) / K;

eqs =@(x) [...
-epsilon^2 + 1/x(2)^2 * (2 / (k+1) * (1 + (k-1) / 2 * x(2)^2 ))^((k+1) / (k-1));...
-F + x(1) * Ve + x(3) * epsilon * (pe - pa);...
-x(1) + pc * x(3) / Cstar
];

x = fsolve(eqs, [0.4838; 2.9355; 9.2e-3], options);

mp = x(1)
Me = x(2)
At = x(3);
dt = A2diam(At)
Ae = epsilon * At;
de = A2diam(Ae)
