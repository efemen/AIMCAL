clear; clc;

% load TrueModelData.mat
load TrueModelDataMagm.mat

% save_path = "./plots/case1_p1_igrf/";
uf = UtilityFunctions();
% s = rng(501); % Setting random number generator seed for repeatability.

%% Sensor reading model

N = length(eulers_deg);

B_tam = zeros(N, 3);          % Model magnetometer measurements.
dcms = eul2rotm(eulers_deg);  % ECI to Body.

n_tam = 300;                  % noise, nT
b_true = [5000; 3000; 6000];  % Constant bias, nT


D_true = [0.05, 0.05, 0.05;   % True D Matrix.
          0.05, 0.1, 0.05;
          0.05, 0.05, 0.05];

x_true = [b_true; D_true(1,1); D_true(2,2); D_true(2,2); D_true(1, 2); D_true(1, 3); D_true(2, 3)];

for i = 1:N
    B_tam(i, :) = ((eye(3) + D_true) \ (dcms(:, :, i) * B_true(i, :)' + b_true + n_tam * randn(3, 1)))';
end


%% EKF Parameters

x = zeros(N, 9);    % EKF estimation vector.
z = zeros(N, 1);    % Measurements
P = zeros(9, 9, N); % Covariances

% P_k =  [1000 * eye(3), zeros(3, 6);
%     zeros(6,3), 0.001 * eye(6)];  

P_k =  [3000 * eye(3), zeros(3, 6);
        zeros(6,3), 1 * eye(6)];  

R = n_tam^2;

b0 = 2.8 * b_true;
D0 = 2.8 * D_true;



D_k0 = [D0(1,1), D0(2,2), D0(2,2), D0(1, 2), D0(1, 3), D0(2, 3)]';

x_k = [b0; D_k0];

x(1, :) = x_k';


%% EKF Execution

for i = 1:N-1
    % Parameters
    P(:,:,i) = P_k;
    b = x_k(1:3);                               % Obtain the bias vector from the current state.
    D = diag(x_k(4:6)) + squareform(x_k(7:9));  % Obtain the scale factor matrix from the current state.

    B = B_tam(i + 1, :)';

    E = 2*D + D^2;
    S = [B(1)^2, B(2)^2, B(3)^2, 2*B(1)*B(2), 2*B(1)*B(3), 2*B(2)*B(3)];
    c = (eye(3) + D) * b;
    
    dEdD = [2*(1 + D(1,1)), 0, 0, 2 * D(1,2), 2 * D(1,3), 0;
            0, 2*(1 + D(2,2)), 0, 2*D(1,2), 0 , 2*D(2,3);
            0, 0, 2*(1 + D(3,3)), 0, 2 * D(1,3), 2 * D(2,3);
            D(1,2), D(1,2), 0, 2 + D(1,1) + D(2,2), D(2,3), D(1,3);
            D(1,3), 0, D(1,3), D(2,3),  2 + D(1,1) + D(3,3), D(1,2); 
            0, D(2,3), D(2,3), D(1,3), D(1,2),  2 + D(2,2) + D(3,3)];

    J = [B(1) * b(1), B(2) * b(2), B(3) * b(3), ...
        B(1)*b(2) + B(2)*b(1), B(1)*b(3) + B(3)*b(1), B(2)*b(3) + B(3)*b(2)];

    H = [2*B' * (eye(3) + D) - 2*b', -S * dEdD + 2 * J];
    
    % Update
    z_k = norm(B_tam(i, :))^2 - norm(B_true(i, :))^2;     % Observation.
   
    K_k = (P_k * H') / (H * P_k * H' + R);                % Kalman gain update.
    x_k = x_k + K_k * (z_k - H * x_k);                    % State update.
    P_k = (eye(9) - K_k * H) * P_k;                       % Covariance matrix update. 
    
    x(i + 1, :) = x_k';
end


x_err = -x + x_true';

rmse = sqrt(sum((x - repmat(x_true', N, 1)).^2) / N)


figure(1)
title("Bias Estimate")

subplot(3, 1, 1)
plot(T, x_err(:, 1), "DisplayName", "b_1 (nT)", "Color", "#ff3366", "LineWidth", 2.5)
ylabel("b_1 Error (nT)")
ylim([-abs(max(x_err(:, 1))), abs(max(x_err(:, 1)))])
grid on
hold on

subplot(3, 1, 2)
plot(T, x_err(:, 2), "DisplayName", "\theta (^o)", "Color", "#ff7f11", "LineWidth", 2.5)
ylabel("b_2 Error (nT)")
ylim([-abs(max(x_err(:, 2))), abs(max(x_err(:, 2)))])
grid on
subplot(3, 1, 3)

plot(T, x_err(:, 3), "Color", "#011627", "LineWidth", 2.5)
ylabel("b_3 Error (nT)")
ylim([-abs(max(x_err(:, 3))), abs(max(x_err(:, 3)))])
xlabel("Time (s)")
grid on
fontsize(15, "points")

set(gcf,'position',[0,0, 1280, 750])
exportgraphics(gcf, save_path + "b_error.png");

figure(2)
title("D Matrix")
subplot(3, 1, 1)
plot(T, x_err(:, 4), "DisplayName", "\omega_1", "Color", "#ffa400", "LineWidth", 2.5)
ylabel("D_{1,1} Error (-)")
ylim([-abs(max(x_err(:, 4))), abs(max(x_err(:, 4)))])
grid on
hold on
subplot(3, 1, 2)
plot(T, x_err(:, 5), "DisplayName", "\omega_2", "Color", "#009ffd", "LineWidth", 2.5)
ylim([-abs(max(x_err(:, 5))), abs(max(x_err(:, 5)))])
ylabel("D_{2,2} Error (-)")
grid on

subplot(3, 1, 3)
plot(T, x_err(:, 6), "DisplayName", "\omega_3", "Color", "#2a2a72", "LineWidth", 2.5)
ylabel("D_{3,3} Error (-)")
grid on
ylim([-abs(max(x_err(:, 6))), abs(max(x_err(:, 6)))])
set(gcf,'position',[0,0, 1280, 750])
xlabel("Time (s)")
fontsize(15, "points")
exportgraphics(gcf, save_path + "d_error.png");


figure(3)
subplot(3, 1, 1)
plot(T, x_err(:, 7), "DisplayName", "\omega_1", "Color", "#ffa400", "LineWidth", 2.5)
ylabel("D_{1,2} Error (-)")
grid on
hold on
ylim([-abs(max(x_err(:, 7))), abs(max(x_err(:, 7)))])


subplot(3, 1, 2)
plot(T, x_err(:, 8), "DisplayName", "\omega_2", "Color", "#009ffd", "LineWidth", 2.5)
ylim([-abs(max(x_err(:, 8))), abs(max(x_err(:, 8)))])
ylabel("D_{1,3} Error (-)")
grid on

subplot(3, 1, 3)
plot(T, x_err(:, 9), "DisplayName", "\omega_3", "Color", "#2a2a72", "LineWidth", 2.5)
ylabel("D_{2,3} Error (-)")
ylim([-abs(max(x_err(:, 9))), abs(max(x_err(:, 9)))])
set(gcf,'position',[0,0, 1280, 750])
fontsize(15, "points")
xlabel("Time (s)")

exportgraphics(gcf, save_path + "d2_error.png");

figure(4)
title("Bias")

subplot(3, 1, 1)
plot(T, x(:, 1) ,"DisplayName", "Estimate", "Color", "r", "LineWidth", 2.5)
hold on
plot(T, ones(N, 1) * b_true(1),"DisplayName", "True", "Color", "k", "LineWidth", 2.5)
ylabel("b_1 (nT)")
ylim([-abs(max(x_err(:, 1))), abs(max(x_err(:, 1)))])
grid on
hold on
legend("Location", "best")

subplot(3, 1, 2)
plot(T, x(:, 2), "DisplayName", "Estimate", "Color", "r", "LineWidth", 2.5)
hold on
plot(T, ones(N, 1) * b_true(2), "DisplayName", "True", "Color", "k", "LineWidth", 2.5)
ylabel("b_2 (nT)")
ylim([-abs(max(x_err(:, 2))), abs(max(x_err(:, 2)))])
grid on
legend("Location", "best")

subplot(3, 1, 3)

plot(T, x(:, 3), "DisplayName", "Estimate", "Color", "r", "LineWidth", 2.5)
hold on
plot(T, ones(N, 1) * b_true(3), "DisplayName", "True", "Color", "k", "LineWidth", 2.5)

ylabel("b_3 (nT)")
ylim([-abs(max(x_err(:, 3))), abs(max(x_err(:, 3)))])
xlabel("Time (s)")
grid on
fontsize(15, "points")
legend("Location", "best")
set(gcf,'position',[0,0, 1280, 750])
exportgraphics(gcf, save_path + "b.png");

figure(5)
title("D Matrix")
subplot(3, 1, 1)
plot(T, x(:, 4), "DisplayName", "Estimate", "Color", "r", "LineWidth", 2.5)
hold on
plot(T, ones(N, 1) * D_k0(1), "DisplayName", "True", "Color", "k", "LineWidth", 2.5)


ylabel("D_{1,1} (-)")
ylim([-abs(max(x_err(:, 4))), abs(max(x_err(:, 4)))])
grid on
legend("Location", "best")

subplot(3, 1, 2)
plot(T, x(:, 5), "DisplayName", "Estimated", "Color", "r", "LineWidth", 2.5)
hold on
plot(T,  ones(N, 1) *D_k0(2), "DisplayName", "True", "Color", "k", "LineWidth", 2.5)
ylim([-abs(max(x_err(:, 5))), abs(max(x_err(:, 5)))])
ylabel("D_{2,2} (-)")
grid on
legend("Location", "best")

subplot(3, 1, 3)
plot(T, x(:, 6), "DisplayName", "Estimated", "Color", "r", "LineWidth", 2.5)
hold on
plot(T, ones(N, 1) * D_k0(3), "DisplayName", "True", "Color", "k", "LineWidth", 2.5)
ylabel("D_{3,3} (-)")
grid on
ylim([-abs(max(x_err(:, 6))), abs(max(x_err(:, 6)))])
set(gcf,'position',[0,0, 1280, 750])
fontsize(15, "points")
legend("Location", "best")
exportgraphics(gcf, save_path + "d.png");


figure(6)
subplot(3, 1, 1)
plot(T, x(:, 7), "DisplayName", "Estimated", "Color", "r", "LineWidth", 2.5)
hold on
plot(T,  ones(N, 1) * D_k0(4), "DisplayName", "True", "Color", "k", "LineWidth", 2.5)
ylabel("D_{1,2} (-)")
grid on
hold on
ylim([-abs(max(x_err(:, 7))), abs(max(x_err(:, 7)))])
legend("Location", "best")

subplot(3, 1, 2)
plot(T, x(:, 8), "DisplayName", "Estimated", "Color", "r", "LineWidth", 2.5)
hold on
plot(T,  ones(N, 1) * D_k0(5), "DisplayName", "True", "Color", "k", "LineWidth", 2.5)
ylim([-abs(max(x_err(:, 8))), abs(max(x_err(:, 8)))])
ylabel("D_{1,3} (-)")
grid on
legend("Location", "best")

subplot(3, 1, 3)
plot(T, x(:, 9), "DisplayName", "Estimated", "Color", "r", "LineWidth", 2.5)
hold on
plot(T,  ones(N, 1) * D_k0(6), "DisplayName", "True", "Color", "k", "LineWidth", 2.5)
ylabel("D_{2,3} (-)")
ylim([-abs(max(x_err(:, 9))), abs(max(x_err(:, 9)))])
set(gcf,'position',[0,0, 1280, 750])
fontsize(15, "points")
legend("Location", "best")
grid on
xlabel("Time (s)")
exportgraphics(gcf, save_path + "d2.png");


