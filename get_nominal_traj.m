function [xnom, ynom] = get_nominal_traj(xinit, tvec, u0, L)    
    [t, xnom] = ode45(@(t,xnom) simulate_nominal_nonlinear_states(t, xnom, u0, L), tvec, xinit);
    
    ynom = simulate_nominal_nonlinear_measurements(xnom);
    xnom = xnom'; ynom = ynom';
end

function dxdt = simulate_nominal_nonlinear_states(t, x, u, L)
    dxdt = [u(1)*cos(x(3)); 
            u(1)*sin(x(3)); 
            u(1)*tan(u(2))/L; 
            u(3)*cos(x(6)); 
            u(3)*sin(x(6)); 
            u(4)];
end

function y = simulate_nominal_nonlinear_measurements(x)
    y = [(atan2((x(:,5)-x(:,2)), (x(:,4)-x(:,1)))-x(:,3)), ...
         sqrt((x(:,1)-x(:,4)).^2+(x(:,2)-x(:,5)).^2), ...
         (atan2((x(:,2)-x(:,5)), (x(:,1)-x(:,4)))-x(:,6)), ...
         x(:,4), ...
         x(:,5)];
end