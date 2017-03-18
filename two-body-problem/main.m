%% FILE main.m

clear;
clc;
close;

%% ODE

% Initial conditions

%    [r1x v1x r1y v1y r1z v1z r2x v2x r2y v2y r2z v2z]'
CI = [0   0   0   0   0   1   20  0   0   4   0   0  ]'; 

% ODE Solver

[t, r] = ode45(@bodiesfun, [0 50], CI);

% Barycenter

r1 = [r(:, 1) r(:, 3) r(:, 5)];
r2 = [r(:, 7) r(:, 9) r(:, 11)];
rb = (r2  - r1) * 1/21 + r1; % r21*m2/(m1+m2) + r1

% Plot

plot3(r(:, 1), r(:, 3), r(:, 5), 'r--');
hold on
plot3(CI(1), CI(3), CI(5), 'r*');
plot3(r(:, 7), r(:, 9), r(:,11), 'b--');
plot3(CI(7), CI(9), CI(11), 'bo');
plot3(rb(:, 1), rb(:, 2), rb(:, 3), 'g-');
