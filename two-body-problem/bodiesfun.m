%% FILE bodiesfun.m

function rp = bodiesfun(t, r) % r is the position

rp = zeros(12, 1); % Output vertical vector

m1 = 20; % kg
m2 = 1; % kg
G = 100; % Gravity constant

% From http://www.math.purdue.edu/~walther/teach/MA366labs/ode45.pdf
% We express the system of multiple order differential equations
% with first order equations/simple derivatives

r12 = [r(1); r(3); r(5)] - [r(7); r(9); r(11)];

F12 = -G*r12/norm(r12)^3; % N/kg2

a1x = F12(1)*m2; % Simplified expression
a1y = F12(2)*m2;
a1z = F12(3)*m2;
a2x = -F12(1)*m1;
a2y = -F12(2)*m1;
a2z = -F12(3)*m1;

% Index of r references
% Odd index of r references positions
% Even index of r references velocities
% rp is the derivative of r of the magnitude the index specifies

%% FIRST BODY (1)
% For the x axis

rp(1) = r(2); % Derivative of the position is the velocity
rp(2) = a1x; % Derivative of velocity is acceleration

% Repeat for y, z axes

rp(3) = r(4);
rp(4) = a1y;

rp(5) = r(6);
rp(6) = a1z;

%% SECOND BODY (2)
% Same procedure

rp(7) = r(8);
rp(8) = a2x;

rp(9) = r(10);
rp(10) = a2y;

rp(11) = r(12);
rp(12) = a2z;

end
