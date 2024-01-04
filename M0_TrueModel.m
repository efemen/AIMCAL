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

X_SC = zeros(N, 3);                  % Position matrix, ECI, km
X_ecef = X_SC;                       % Position matrix, ECEF, km
V_SC = X_SC;                         % Velocity matrix, km / s
LLA_SC = X_SC;                       % Lat-lon-alt matrix, meters
eulers = X_SC;                       % Euler angles, yaw, pitch, roll, rad
w = X_SC;                            % Body rates, rad / s

B_true = X_SC;                            % Magnetic field matrix, ECI frame, nT

X_SC(1, :) = pos;                    % Initial position, km
V_SC(1, :) = vel;                    % Initial velocity, km / s
eulers(1, :) = [0, 0, 0];            % Initial attitude, yaw, pitch, roll, rad
w(1, :) = deg2rad([-0.1, 0.5, 0.5]); % Initial body rates, rad / s


%% EOM

a = @(X) -earth.mu * X / norm(X)^3; % Acceleration expression
w_dot_fun = @(w) [(-(I3 - I2) * w(2) * w(3)) / I1; 
                  (-(I1 - I3) * w(3) * w(1)) / I2; 
                  (-(I2 - I1) * w(1) * w(2)) / I3]'; % Torque free rotational dynamics.

%% Solution loop

for i = 1:N-1
    t = T(i);
    sdt = sdt_vector(i);

    [X_SC, V_SC] = uf.RK4(a, dt, X_SC, V_SC, i); % Propagate orbit
    [eulers, w] = F6_RK4(w_dot_fun, dt, eulers, w, i); % Propagate attitude
    
    X_ecef(i, :) = uf.ECI2ECEF(X_SC(i, :), sdt);  % ECI to ECEF conversion.

    LLA_SC(i, :) = ecef2lla(1e3 * X_ecef(i, :), 'WGS84');    % Lat, lon, alt from ECEF
    % B_NED = igrfmagm(LLA_SC(i, 3), LLA_SC(i, 2), LLA_SC(i, 1), decyear(date_vector(i)) - 20); % 2043 is invalid.
    B_NED = wrldmagm(LLA_SC(i, 3), LLA_SC(i, 2), LLA_SC(i, 1), decyear(date_vector(i)) - 20)'; % 2043 is invalid.
    B_true(i, :) = uf.ECEF2ECI(uf.NED2ECEF(B_NED', LLA_SC(i, 1), LLA_SC(i, 2))', sdt);
end

eulers_deg = uf.AngleVectorNormalizerDeg(rad2deg(eulers));
rates_deg = rad2deg(w);

figure(1)
title("Euler Angles")
subplot(3, 1, 1)
plot(T, eulers_deg(:, 1), "DisplayName", "\psi (^o)", "Color", "#ff3366", "LineWidth", 2.5)
ylabel("\psi (^o)")
grid on
hold on
subplot(3, 1, 2)
plot(T, eulers_deg(:, 2), "DisplayName", "\theta (^o)", "Color", "#ff7f11", "LineWidth", 2.5)
ylabel("\theta (^o)")
grid on
subplot(3, 1, 3)
plot(T, eulers_deg(:, 3), "Color", "#011627", "LineWidth", 2.5)
ylabel("\phi (^o)")

xlabel("Time (s)")
grid on
fontsize(15, "points")

set(gcf,'position',[0,0, 1280, 750])

figure(2)
title("Body Rates")
subplot(3, 1, 1)
plot(T, rates_deg(:, 1), "DisplayName", "\omega_1", "Color", "#ffa400", "LineWidth", 2.5)
ylabel("\omega_1 (^o/s)")
grid on
hold on
subplot(3, 1, 2)
plot(T, rates_deg(:, 2), "DisplayName", "\omega_2", "Color", "#009ffd", "LineWidth", 2.5)
ylabel("\omega_2 (^o/s)")
grid on

subplot(3, 1, 3)
plot(T, rates_deg(:, 3), "DisplayName", "\omega_3", "Color", "#2a2a72", "LineWidth", 2.5)
ylabel("\omega_3 (^o/s)")

set(gcf,'position',[0,0, 1280, 750])
fontsize(15, "points")


xlabel("Time (s)")
grid on

figure(3)
title("B_{true}")
subplot(3, 1, 1)
plot(T(1:end-1), B_true(1:end-1, 1), "DisplayName", "\omega_1", "Color", "#ffa400", "LineWidth", 2.5)
ylabel("B_x")
grid on
hold on
subplot(3, 1, 2)
plot(T(1:end-1), B_true(1:end-1, 2), "DisplayName", "\omega_2", "Color", "#009ffd", "LineWidth", 2.5)
ylabel("B_y")
grid on

subplot(3, 1, 3)
plot(T(1:end-1), B_true(1:end-1, 3), "DisplayName", "\omega_3", "Color", "#2a2a72", "LineWidth", 2.5)
ylabel("B_z")

set(gcf,'position',[0,0, 1280, 750])
fontsize(15, "points")


xlabel("Time (s)")
grid on


% save TrueModelData.mat T X_SC B_true eulers_deg
% save TrueModelDataMagm.mat