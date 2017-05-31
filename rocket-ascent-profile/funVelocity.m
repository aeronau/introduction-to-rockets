% FILE funVelocity.m

function rp = funVelocity(t, r, stageInfo)

    % r = [h v k A m phi lambda]

    T = stageInfo.T; % N
    mdot = stageInfo.mdot; % At SL, kg/s
    CD = stageInfo.CD; % From 0.75 to 2.4
    Area = stageInfo.Area; % m2
    Re = 6373e3; % m
    M = 5.9722e24; % kg
    omega = 2*pi/3600/24; % Earth's angular speed, rad/s
    L = 0; % Lift is assumed null, N
    tturn = 5; %time until gravity turn, s

    % CONNECTIVITY TABLE
    % ------------------
    % 1 for h
    % 2 for v
    % 3 for k
    % 4 for A
    % 5 for m
    % 6 for phi
    % 7 for lambda
    % ------------------

    h = r(1);
    v = r(2);
    k = r(3); % gamma is k
    A = r(4);
    m = r(5);
    phi = r(6);
    lambda = r(7);

    % Equations

    dmdt = -mdot;

    dist_to_orig = Re + h;

    g = gravity(M, h, Re);

    [~, ~, density] = exponentialAtmosphere(h);

    D = 1/2 * density * v^2 * Area * CD;
    
    dvdt = T/m - D/m - g*sin(k) + omega^2 * dist_to_orig * (cos(phi))^2 * (sin(k) - cos(k) * tan(phi) * cos(A));
    
    if (t < tturn)
    
        dkdt = 0;
    
        dAdt = 0;
            
    else

        dkdt = 1/v * (L/m - g*cos(k) + v^2/dist_to_orig*cos(k) + 2*omega*v*cos(phi)*sin(A) + omega^2 * dist_to_orig *(cos(phi))^2*(cos(k) + sin(k)*tan(phi)*cos(A)));

    end
    
    dAdt = 0;

    dmdt = -mdot;

    vh = v*sin(k);

    dlambdadt = v*sin(A)*cos(k) / dist_to_orig / cos(phi);

    dphidt = v*cos(A)*cos(k)/dist_to_orig;

    % OUTPUT SETUP
    
    rp(1) = vh;
    rp(2) = dvdt;
    rp(3) = dkdt;
    rp(4) = dAdt;
    rp(5) = dmdt;
    rp(6) = dphidt;
    rp(7) = dlambdadt;

    rp = [vh dvdt dkdt dAdt dmdt dphidt dlambdadt]';

end
