function e_angles_dot = F3_EulerDot(w, e_angles)
% Returns the time rate of change of euler angles w.r.t. instant euler
% angles and angular velocity vector
    yaw = e_angles(1);
    pitch = e_angles(2);
    roll = e_angles(3);
    
    e1 = [0 sin(roll) cos(roll)];
    e2 = [0 cos(roll)*cos(pitch) -sin(roll)*cos(pitch)];
    e3 = [cos(pitch) sin(roll)*sin(pitch) cos(roll)*sin(pitch)];

    e_angles_dot = (1/cos(pitch)) * [e1; e2; e3] * w';
    e_angles_dot = e_angles_dot';
end

