clear; clc; close all
%% Initial state vector of parking orbit
load parking_data.mat

uf = UtilityFunctions();
jdt = juliandate(parking_orbit_date);

earth = CelestialObject("Earth", 5.97217e24, 6371.0084, 1.49598e8, 23.43928, jdt); 
                        % name, mass, planetary radius, heliocentric radius, jdt

elements = OrbitalElements(pos, vel, earth.mu); % Orbital elements

m = 2.5197E+04; % kg, most of the mass belongs to the Centaur Upper Stage.
I = m * [1.2 0 0; 0 1.1 0; 0 0 0.8]; % kg * m^2

I1 = I(1, 1);
I2 = I(2, 2);
I3 = I(3, 3);

%% Discretization

dt = 0.1;                            % seconds
T = 0:dt:elements.period;            % Time vector

date_vector = parking_orbit_date + seconds(T);
jdt_vector = juliandate(date_vector);
sdt_vector = siderealTime(jdt_vector);

N = length(T);                       % Iteration length

X_SC = zeros(N, 3);                  % Position matrix, km
V_SC = X_SC;                         % Velocity matrix, km / s
LLA_SC = X_SC;                       % Lat-lon-alt matrix, meters
eulers = X_SC;                       % Euler angles, yaw, pitch, roll, rad
w = X_SC;                            % Body rates, rad / s

B = X_SC;                            % Magnetic field matrix, ECI frame, nT

X_SC(1, :) = pos;                    % Initial position, km
V_SC(1, :) = vel;                    % Initial velocity, km / s
eulers(1, :) = [0, 0, 0];            % Initial attitude, rad
w(1, :) = deg2rad([-0.1, 0.5, 0.5]); % Initial body rates, rad / s


%% EOMs

a = @(X) -earth.mu * X / norm(X)^3; % Acceleration expression
w_dot_fun = @(w) [(-(I3 - I2) * w(1,2) * w(1,3)) / I1; 
                  (-(I1 - I3) * w(1,3) * w(1,1)) / I2; 
                  (-(I2 - I1) * w(1,1) * w(1,2)) / I3]'; % Torque free rotational dynamics.

for i = 1:N-1
    t = T(i);
    sdt = sdt_vector(i);

    [X_SC, V_SC] = uf.RK4(a, dt, X_SC, V_SC, i); % Propagate orbit
    [eulers, w] = uf.RK4_euler(w_dot_fun, dt, eulers, w, i); % Propagete attitude
    
    x_ecef = uf.ECI2ECEF(X_SC(i, :), sdt);

    LLA_SC(i, :) = ecef2lla(1e3 * x_ecef', 'WGS84');
    B_NED = igrfmagm(LLA_SC(i, 3), LLA_SC(i, 2), LLA_SC(i, 1), decyear(date_vector(i)) - 20); % 2043 is invalid.
    B(i, :) = uf.ECEF2ECI(uf.NED2ECEF(B_NED', LLA_SC(i, 1), LLA_SC(i, 2))', sdt);
end

eulers_deg = uf.AngleVectorNormalizerDeg(rad2deg(eulers));

