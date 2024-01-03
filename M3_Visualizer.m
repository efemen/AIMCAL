clear; clc; close all;

%% Assumptions

% All planetary orbits are coplanar.
% Earth is a sphere.
% No plane change except from launch to equatorial orbit.


%% User Parameters

load parking_data.mat
load TrueModelData.mat

dcms = eul2rotm(eulers_deg);  % ECI to Body.
jdt = juliandate(parking_orbit_date);
sdt = siderealTime(jdt); % Earth map update

%% Object Initilization

uf = UtilityFunctions();

earth = CelestialObject("Earth", 5.97217e24, 6371.0084, 1.49598e8, 23.43928, jdt); % name, mass, planetary radius, heliocentric radius, jdt

%% Setup Geometry and Plots
% Earth Sphere
[X,Y,Z] = sphere(100);

X_E = X * earth.r;
Y_E = Y * earth.r;
Z_E = Z * earth.r;
c_RotX = mean(mean(X_E));
c_RotY = mean(mean(Y_E));
c_Rot = [c_RotX c_RotY 0];

% Ecliptic Plane
n_ecliptic = uf.t_xX(90-earth.tilt, 270) * [0; 0; 1]; % RA = 270 deg, DEC = 66.56 deg

[X_ecliptic, Y_ecliptic] = meshgrid(-8000:2000:8000); % Generate x and y data
Z_ecliptic = -1/n_ecliptic(3) * (n_ecliptic(1)*X_ecliptic + n_ecliptic(2)*Y_ecliptic); % Solve for z data

% Sun Direction
n_sun =  uf.ICRF2ECI((uf.hat(-earth.heliocentric_pos)'));

% Earth Velocity Direction
n_earth_velocity = uf.rodrigues_rot(n_sun, n_ecliptic, -90);

% Plots

figure(1);
earth_map = surf(X_E,Y_E,-Z_E);
earthMap = imread("world_Map.jpg");
set(earth_map,'CData', earthMap,'FaceColor','texturemap',"EdgeColor","none")
hold on
colormap white
axis equal
rotate(earth_map, [0 0 1], sdt)
camva(3)
uf.draw_space();

ecliptic_plane = surf(X_ecliptic, Y_ecliptic, Z_ecliptic);
ecliptic_plane.FaceAlpha = 0.2;
ecliptic_plane.EdgeColor = "white";
ecliptic_plane.EdgeAlpha = 0.2;
hold on

att = quiver3(X_SC(1, 1), X_SC(1, 2), X_SC(1, 3), 1, 0, 0, 1000, "filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "red");
b_field = quiver3(X_SC(1, 1), X_SC(1, 2), X_SC(1, 3), 1, 0, 0, 5E-2, "filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "green");
telemetry = text(X_SC(1, 1), X_SC(1, 2), X_SC(1, 3), "", "Color", "white", "FontSize", 20);

set(gcf, 'Position',  [0, 0, 1920, 1080])

clear X_E Y_E Z_E X Y Z X_ecliptic Y_ecliptic Z_ecliptic



%% Setup Spacecraft Initial Conditions

dt = 0.1; % s
N = length(X_SC);    % Iteration length
T = 0:dt:N*dt;     % Time matrix


earth_w = 7.292115e-5; % rad / s

for i = 1:6000:N
    T(i + 1) = T(i) + dt;
    [ra_r, dec_r] = uf.ECI2raDec(X_SC(i, :));

    current_att = dcms(:, :, i) * [1; 0; 0];
    att.XData = X_SC(i, 1);
    att.YData = X_SC(i, 2);
    att.ZData = X_SC(i, 3);
    att.UData = current_att(1);
    att.VData = current_att(2);
    att.WData = current_att(3);

    b_field.XData = X_SC(i, 1);
    b_field.YData = X_SC(i, 2);
    b_field.ZData = X_SC(i, 3);
    b_field.UData = B_true(i, 1);
    b_field.VData = B_true(i, 2);
    b_field.WData = B_true(i, 3);

   % Plot current position.
    figure(1)
    plot3(X_SC(i,1), X_SC(i, 2), X_SC(i, 3),"Color","#FF3131");
    rotate(earth_map, [0 0 1], rad2deg(earth_w*dt), c_Rot)
    view(ra_r + sdt + 110, dec_r+5);
    telemetry.String = "T+ " + string(i * dt - 0.1) + " s";
    telemetry.Position = X_SC(i ,:) * 1.5;

    pause(2)
    exportgraphics(gca, "./plots/" + string(i) + ".png")

end

