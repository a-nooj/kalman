function [x, y,dx] = simulate_linearized_system(xnom, ynom, x0, L, u, delta_T, num_timesteps, dx0)
    dx = [dx0];
    x = [x0+dx0];
    y = [NaN; NaN; NaN; NaN; NaN];

    for k=2:num_timesteps+1
        [F, G, ~, ~] = get_discrete_mats(xnom(:, k-1), u, L, delta_T);
        
        dx_next = F*dx(:, k-1);
        dx = [dx, dx_next];
        x_next = xnom(:,k) + dx_next;
        x = [x, x_next];
        
        [~, ~, H, ~] = get_discrete_mats(xnom(:,k), u, L, delta_T);
        ythis = ynom(:, k-1) + H*dx(:, k);
        y = [y, ythis];
    end
end