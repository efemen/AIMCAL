classdef UtilityFunctions
    % This file organizes the functions that are utilized in this work in a
    % document.
    
    methods
        function obj = UtilityFunctions()
            disp("UtilityFunctions initialized.")
        end
        
        function e_angles_normalized = AngleVectorNormalizerDeg(obj, e_angles)
        % Loops through the whole vector to normalize angles.
            n = size(e_angles);
            e_angles_normalized = zeros(n(1), n(2));
            
            for i = 1:n(1)
                for j = 1:n(2)
                    e_angles_normalized(i, j) = obj.AngleNormalizerDeg(e_angles(i, j));
                end
            end
        end
        
        function normalizedAngle = AngleNormalizerDeg(obj, angle)
        % Normalizes an angle to a range of [-180, 180]
            if angle <= 180 && angle >= -180
                normalizedAngle = angle;
            end
   
            if angle > 180
                int_multiples = floor(angle / 360);
                normalizedAngle = angle - 360 * int_multiples;
                if normalizedAngle > 180
                    normalizedAngle = normalizedAngle - 360;
                end
            elseif angle < -180
                int_multiples = floor(angle / 360);
                normalizedAngle = angle - 360 * int_multiples;
                if normalizedAngle > 180
                    normalizedAngle = normalizedAngle - 360;
                end
            end
        end


        function e_angles_dot = euler_dot(obj, w, e_angles)
        % Returns the time rate of change of euler angles w.r.t. instant euler
        % angles and angular velocity vector
            yaw = e_angles(1);
            pitch = e_angles(2);
            roll = e_angles(3);
            
            e1 = [0 sin(roll) cos(roll)];
            e2 = [0 cos(roll)*cos(pitch) -sin(roll)*cosd(pitch)];
            e3 = [cos(pitch) sin(roll)*sin(pitch) cos(roll)*sind(pitch)];
        
            e_angles_dot = (1/cos(pitch)) * [e1; e2; e3] * w';
            e_angles_dot = e_angles_dot';
        end

        function x_ecef = ECI2ECEF(obj, X, sdt)
            x_ecef = [cosd(sdt) sind(sdt) 0;
                      -sind(sdt) cosd(sdt) 0;
                      0,    0,      1] * X';
        end

        function x_ecef = NED2ECEF(obj, x_ned, lat, lon)
            dcm =  [-sind(lat) * cosd(lon), -sind(lat) * sind(lon),  cosd(lat);...
                   -sind(lon),            cosd(lon),        0;...
                   -cosd(lat) * cosd(lon), -cosd(lat) * sind(lon), -sind(lat)]';
            x_ecef = dcm * x_ned;
        end

        function x_eci = ECEF2ECI(obj, X, sdt)
            x_eci = [cosd(sdt) sind(sdt) 0;
                      -sind(sdt) cosd(sdt) 0;
                      0,    0,      1]' * X';
        end


        function v_hat = hat(obj, v)
            v_hat = v / norm(v);
            if isnan(v_hat)
                v_hat = zeros(size(v));
            end
        end

        function angle = angle_between(obj, v1, v2)
            angle = acosd(dot(v1, v2) / (norm(v1) * norm(v2)));
            if isnan(angle)
                angle = 0;
            end
        end
        
        function [X_next, V_next] = RK4(obj, dydx, dt, X_SC, V_SC, i)
        % RK4 numerical solver
            dv1 = dt * dydx(X_SC(i,:));
            dx1 = dt * V_SC(i,:);
        
            dv2 = dt * dydx(X_SC(i,:) + 0.5 * dx1);
            dx2 = dt * (V_SC(i,:) + 0.5 * dv1);
        
            dv3 = dt * dydx(X_SC(i,:) + 0.5  * dx2);
            dx3 = dt * (V_SC(i,:) +  0.5 * dv2);
        
            dv4 = dt * dydx(X_SC(i,:) + dx3);
            dx4 = dt * (V_SC(i,:) + dv3);
        
            V_SC(i + 1,:) = V_SC(i,:) +  (dv1 + 2*dv2 + 2*dv3 + dv4) / 6;
            X_SC(i + 1,:) = X_SC(i,:) + (dx1 + 2*dx2 + 2*dx3 + dx4) / 6;
        
            X_next = X_SC;
            V_next = V_SC;
        end


        function out = ResultPlot(obj, T, x, P, true_val)
        N = length(T);
        var_names = ["b_1 (nT)", "b_2 (nT)", "b_3 (nT)","D_{1,1}  (-)", ...
                    "D_{2,2} (-)", "D_{3,3}(-)", "D_{1,2}(-)", "D_{1,3}(-)", ...
                    "D{2,3}(-)"];
    
            for i = 1:9
                figure(i)
                set(gcf, "Position", [100 100 1000 600])
                
                plot(T, x(:, i) - true_val(i), "Color", "k");
                hold on
                grid on
                xlabel("Time (s)")
                ylabel(var_names(i))
        
                bounds = zeros(N, 1);
                for j = 1:N
                    bounds(j)  = 3 * sqrt(P(i, i, j));
                end
                
                % ylim([min(-bounds) max(bounds)])
                % plot(T(1:end-1), bounds(1:end-1), "--", "Color", "red")
                % plot(T(1:end-1), -bounds(1:end-1), "--", "Color", "red")
                % xlim([0 100])
            end

        end
           function xX = t_xX(obj, lat, sid)
            % Topocentric to geocentric equatorial
            xX = zeros(3,3);
            xX(1,:) = [-sind(sid) -sind(lat)*cosd(sid) cosd(lat)*cosd(sid)];
            xX(2,:) = [cosd(sid) -sind(lat)*sind(sid) cosd(lat)*sind(sid)];
            xX(3,:) = [0 cosd(lat) sind(lat)];
        end

        function [ra, dec] = ECI2raDec(obj, R)
            % Convert ECI coordinates to right ascension and declination.
            R_hat = R/norm(R);
            dec = asind(R_hat(3));
            if R(2) > 0
                ra = acosd(R_hat(1) / cosd(dec));
            else
                ra = 360 - acosd((R_hat(1) / cosd(dec)));
            end
        end

        function X_ECI = ICRF2ECI(obj, X)
            X_ECI = [1, 0, 0;...
          0, cosd(23.44), -sind(23.44);...
          0, sind(23.44),  cosd(23.44)] * X;
        end

        function X_ICRF = ECI2ICRF(obj, X)
          X_ICRF = [1, 0, 0;...
          0, cosd(-23.44), -sind(-23.44);...
          0, sind(-23.44),  cosd(-23.44)] * X;
        end


        function v_rot = rodrigues_rot(obj, v, k, angle)
            v_rot = v * cosd(angle) + cross(k,v) * sind(angle) + k * dot(k, v) * (1 - cosd(angle));
        end


        function dark = draw_space(obj)
            set(gcf,'Color','black');
            set(gca,'Color','black');
            set(gca, 'GridColor', 'white'); 
            set(gca, 'GridAlpha', 0.5);
            set(gca, "XColor", "white");
            set(gca, "YColor", "white");
            set(gca, "ZColor", "white");
            dark = 1;
        end

    end
end

