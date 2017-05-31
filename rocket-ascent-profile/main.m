% FILE main.m

close;
clear;
clc;

%% KNOWN INITIAL VALUES
% Magnitude = Value; % Equivalence in array

initial_position = 0; % h
initial_velocity = 0; % v
k = deg2rad(89.5); % k
azimuth = deg2rad(90); % A
latitude = deg2rad(5.232); % phi
longitude = deg2rad(-52.776); % lambda


%% ROCKET STAGE INFORMATION
% Set integration limits

numSequences = 5; % Including coasting
         % [time_start time_to_1st_seq time_to_2nd_seq time_to_3rd_seq time_to_4th_seq  time_to_5th_seq];
timeMark = [0          115             202             226             231              354            ]; % Vector of length numSequences+1 or only reference time if duration of each sequence is specified in timeDelta, s
timeDelta = diff(timeMark); % Duration of each rocket sequence, s

% Set known values for each stage
% mdot, T, CD, Area

mdot = [87710/109.9, 23814/77.1, 0, 10567/119.6, 10567/119.6]; % Fuel consumption values for each sequence (0 implies coasting happening). It is assumed constant throughout the burnout time and equal to the propellant_mass/burnout_time. kg/s
T = [2261e3 1196e3 0 225e3 225e3]; % Same for thrust, N
CD = [1 1 1 1 1]; % Same for CD
Area = pi/4*[3 1.9 1.9 1.9 1.9].^2; % Same for reference area, which is the cross section of the biggest stage. Array of diameters, m2

% Individual gross mass for each stage

mFairing = 540; % kg

mOrbital = 688 + 77; % Mass for the fairing, upper stage and the PLA 937 VG adapter. This mass is always present during the ascent. kg

mSeq = [96243 26300 mFairing 12000 mOrbital]; % kg

% IMPORTANT: Coasting happening ONLY from 202 to 226s

m = [mSeq(1) + mSeq(2) + mSeq(3) + mSeq(4) + mSeq(5),... % All stages
     mSeq(2) + mSeq(3) + mSeq(4) + mSeq(5),... % No 1st stage
     mSeq(3) + mSeq(4) + mSeq(5),... % No 1st and 2nd stage
     mSeq(3) + mSeq(4) + mSeq(5),... % No 1st and 2nd stage during coasting
     mSeq(4) + mSeq(5)... % No 1st and 2nd stage and no fairing
    ]; % Mass of rocket at each timemark, Kg

% Initial mass is  the sum of all stages
% http://www.arianespace.com/wp-content/uploads/2015/09/Vega-Users-Manual_Issue-04_April-2014.pdf

initial_mass = sum(mSeq); % Sum non-repeated masses of each sequence, m

for i = 1:numSequences

    % Build an array of structures containing info for each sequence

    stages(i) = struct('m', m(i),...
                       'mdot', mdot(i),...
                       'T', T(i),...
                       'CD', CD(i),...
                       'Area', Area(i));

end

%% ODE SOLVER CALL
% init

odeset = []; % Options will not be used
tSave = []; % Save timemarks
rSave = []; % Save values according to timemarks

% INDEX FOR ARRAY CI
% ------------------
% 1 for h
% 2 for v
% 3 for k
% 4 for A
% 5 for m
% 6 for phi
% 7 for lambda
% ------------------

r = [initial_position initial_velocity k azimuth initial_mass latitude longitude]; % init r vector

if m(1) ~= initial_mass % Basic mass checkup
    error('MASS IS NOT CONSISTENT: Check that in m vector, the first element includes the mass of all stages.');
end

for i = 1:numSequences;

    CI = r(end, :);

    CI(5) = m(i); % Replace m at each iter because the mass function is not continuous (mass of each stage is lost at their separation). Not pretty but faster to write than overwriting all variables except m

    [t, r] = ode45(@funVelocity, [timeMark(1) timeDelta(i)], CI, odeset, stages(i));

    % Append results (avoids overriding values)

    tSave = [tSave; t + timeMark(i)]; % Because we are integrating from 0 to the duration of each stage, we have to translate the points to when there is a change of stage
    rSave = [rSave; r];

end

subplot(4, 2, 2);
plot(tSave,rSave(:,2));
xlabel('Time (s)');
ylabel('Velocity (m/s)');

subplot(4, 2, 1);
plot(tSave,rSave(:, 1));
xlabel('Time (s)');
ylabel('Height (m)');

subplot(4, 2, [3, 4]);
plot(tSave,rSave(:,5));
xlabel('Time (s)');
ylabel('Mass (kg)');
tPoints = ismember(tSave, timeMark);
rPoints = rSave(tPoints, 5);
hold on
plot(tSave(tPoints), rPoints, 'g*');
text(timeMark(1), rPoints(1), 'init');
text(timeMark(2), rPoints(2), '1/2 sep.');
text(timeMark(3), rPoints(4), '2/3 sep.', 'HorizontalAlignment', 'right');
text(timeMark(4), rPoints(6), 'Coasting','VerticalAlignment','top', 'HorizontalAlignment', 'center');
text(timeMark(5), rPoints(8), 'Fairing sep.', 'VerticalAlignment','bottom');
text(timeMark(end), rPoints(end), '3/4 sep.','VerticalAlignment','bottom', 'HorizontalAlignment', 'center');

subplot(4, 2, 5);
plot(tSave,rSave(:,7));
xlabel('Time (s)');
ylabel('Longitud (rad)');

subplot(4, 2, 6);
plot(tSave,rSave(:,6));
xlabel('Time (s)');
ylabel('Latitude (rad)');

subplot(4, 2, [7,8]);
plot(tSave,rSave(:,3));
xlabel('Time (s)');
ylabel('Gamma (rad)');
