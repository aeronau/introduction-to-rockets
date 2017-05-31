%% FILE exponentialAtmosphere.m
% Returns T, P, rho at a height h

function [T, P, rho] = exponentialAtmosphere(h) % h is a vector of altitudes at which T, P, rho are desired

    PSL = 1e5; % Pa
    rhoSL = 1.225; % kg/m3

    n = length(h);

    % Temperatures for each altitude range, in K

    T(1:n) = 210;
    T(h<=5e5) = 260;
    T(h<=86e3) = 273;
    T(h<=11e3) = 290;

    % Constants of height for each altitude range, in m

    H(1:n) = 6e3;
    H(h<=5e5) = 7.61e3;
    H(h<=86e3) = 8e3;
    H(h<=11e3) = 8.5e3;

    P = PSL.*exp(-h./H');
    rho = rhoSL.*exp(-h./H');

    T = T';

end
