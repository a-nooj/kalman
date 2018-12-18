function [x, y] = generate_noisy_ground_truth_data(u0, L, q, R, t0, t1, xinit)
    [t, x] = ode45(@(t,x) simulate_noisy_nonlinear_states(t, x, u0, L, q), ...
        [t0 t1], xinit);
    x = x(end,:);

    y = simulate_noisy_nonlinear_measurements(x, R);
    x = x'; y = y';
end

function dxdt = simulate_noisy_nonlinear_states(t, x, u, L, noise_mat)
    dxdt = [u(1)*cos(x(3)); 
            u(1)*sin(x(3)); 
            u(1)*tan(u(2))/L; 
            u(3)*cos(x(6)); 
            u(3)*sin(x(6)); 
            u(4)] + noise_mat;
end

function y = simulate_noisy_nonlinear_measurements(x, R)
    y = [(atan2((x(5)-x(2)), (x(4)-x(1)))-x(3)), ...
         sqrt((x(1)-x(4)).^2+(x(2)-x(5)).^2), ...
         (atan2((x(2)-x(5)), (x(1)-x(4)))-x(6)), ...
         x(4), ...
         x(5)] + mvnrnd(zeros(1, length(R)), R);
end