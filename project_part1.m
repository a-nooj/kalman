function [] = project_part1()
    [delta_T, L, x0, u0, num_timesteps, ~, ~] = get_common_variables();
    [F, G, H, ~] = get_discrete_mats(x0, u0, L, delta_T);
    
    num_states = length(x0); num_inputs = length(u0);
    
    % hw8 2b
    
    observable = rank(obsv(F,H))==num_states;
    stable = eig(F);
    controllable = rank(ctrb(F,G))==num_states;
    
    given = load('cooplocalization_finalproj_KFdata.mat');
    Qtrue = given.Qtrue;
    Rtrue = given.Rtrue;
    yreal = given.ydata;
    tvec = given.tvec;
    measLabels = given.measLabels;
    
    % hw8 2c
    
    % nominal traj
    [xnom, ynom] = get_nominal_traj(x0,tvec, u0, L);
    
    % nonlinear traj
    dx = [0; 1; 0; 0; 0; 0.1];  % initial state perturbation
    [t, xnonlin] = ode45(@(t,xnonlin) simulate_nonlinear_states(t, xnonlin, u0, L, false, 0), 0:delta_T:num_timesteps*delta_T, x0+dx);
    xnonlin(:,3) = wrapToPi(xnonlin(:,3));
    xnonlin(:,6) = wrapToPi(xnonlin(:,6));
    %off by one?
    ynonlin = simulate_nonlinear_measurements(xnonlin, false, 0);
    xnonlin = xnonlin'; ynonlin = ynonlin';
    
    [xlin, ylin,dx] = simulate_linearized_system(xnom, ynom, x0, L, u0, delta_T, num_timesteps, dx);
    xlin(3,:) = wrapToPi(xlin(3,:));
    xlin(6,:) = wrapToPi(xlin(6,:));
    %dx(3,:) = wrapToPi(dx(3,:));
    %dx(6,:) = wrapToPi(dx(6,:));
    ynonlin(:,1)=[nan;nan;nan;nan;nan];
    %plot_states_measurements(t, xnonlin, ynonlin);
    %plot_states_measurements(t, xlin, ylin);
    plot_states_measurements(t,dx,ynonlin-ylin);
end