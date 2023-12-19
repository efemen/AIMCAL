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
    
       function [e_angles_next, w_next] = RK4_euler(obj, dwdt, dt, e_angles, w, i)
            % RK4 numerical solver
        
            dw1 = dt * dwdt(w(i,:));
            d_euler1 = dt * obj.euler_dot(w(i,:), e_angles(i,:));
        
            dw2 = dt * dwdt(w(i,:) + 0.5 * dw1);
            d_euler2 = dt * obj.euler_dot(w(i,:) + 0.5 * dw1, e_angles(i,:) + 0.5 * d_euler1);
        
            dw3 = dt * dwdt(w(i,:) + 0.5  * dw2);
            d_euler3 = dt * obj.euler_dot(w(i,:) + 0.5 * dw2, e_angles(i,:) + 0.5 * d_euler2);
        
            dw4 = dt * dwdt(w(i,:) + dw3);
            d_euler4 = dt * obj.euler_dot(w(i,:) + dw3, e_angles(i,:) +  d_euler3);
        
            w(i + 1,:) = w(i,:) +  (dw1 + 2*dw2 + 2*dw3 + dw4) / 6;
            e_angles(i + 1,:) = e_angles(i,:) + (d_euler1 + 2*d_euler2 + 2*d_euler3 + d_euler4) / 6;
        
            e_angles_next = e_angles;
            w_next = w;
end


    end
end

