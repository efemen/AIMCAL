function [e_angles_next, w_next] = RK4(dwdt, dt, e_angles, w, i)
% RK4 numerical solver

    dw1 = dt * dwdt(w(i,:));
    d_euler1 = dt * F3_EulerDot(w(i,:), e_angles(i,:));

    dw2 = dt * dwdt(w(i,:) + 0.5 * dw1);
    d_euler2 = dt * F3_EulerDot(w(i,:) + 0.5 * dw1, e_angles(i,:) + 0.5 * d_euler1);

    dw3 = dt * dwdt(w(i,:) + 0.5  * dw2);
    d_euler3 = dt * F3_EulerDot(w(i,:) + 0.5 * dw2, e_angles(i,:) + 0.5 * d_euler2);

    dw4 = dt * dwdt(w(i,:) + dw3);
    d_euler4 = dt * F3_EulerDot(w(i,:) + dw3, e_angles(i,:) +  d_euler3);

    w(i + 1,:) = w(i,:) +  (dw1 + 2*dw2 + 2*dw3 + dw4) / 6;
    e_angles(i + 1,:) = e_angles(i,:) + (d_euler1 + 2*d_euler2 + 2*d_euler3 + d_euler4) / 6;

    e_angles_next = e_angles;
    w_next = w;
end