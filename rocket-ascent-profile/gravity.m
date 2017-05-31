% FILE gravity.m
% Outputs val of gravity at an altitude h with respect to the surface of Earth

function g = gravity(M, h, Re) % Mass (kg), Earth radius (m) and altitude h (m) as inputs

    G = 6.67384e-11; % Newton's Universal Gravity Constant, m/s2

    g = G.*M./(Re+h).^2;

end
