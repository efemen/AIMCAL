classdef OrbitalElements
    % Initializing this object with the state vector and the center body
    % returns the orbital elements in properties.
    
    properties
        h_vector
        h
        e_vector
        e
        r_periapsis
        r_apoapsis
        period
        a 
    end
    
    methods
        function obj = OrbitalElements(X_i, V_i, mu)
            % Calculate orbital elements
            obj.h_vector = cross(X_i, V_i);
            obj.h = norm(obj.h_vector);

            V_ri = dot(X_i/norm(X_i), V_i);
            obj.e_vector = (1/mu) * ((norm(V_i)^2 - mu / norm(X_i)) * X_i - norm(X_i) * V_ri * V_i);
            obj.e = norm(obj.e_vector);
            obj.r_apoapsis = (obj.h^2 / mu) ./ (1 + obj.e * cosd(180));
            obj.r_periapsis = (obj.h^2 / mu) ./ (1 + obj.e * cosd(0));
            obj.a = 0.5 * (obj.r_periapsis + obj.r_apoapsis);
            obj.period =  sqrt(obj.a^3 * 4 * pi^2 / mu);
        end
        
    end
end

