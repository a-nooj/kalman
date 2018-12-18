function [] = plot_state_errors(t, state_est_errs, pos_2sig, neg_2sig,sigs)
    figure;

    subplot(6,1,1);
    plot(t, state_est_errs(1,:), 'b'); hold on;
    plot(t, pos_2sig(1,:), 'r--'); hold on;
    plot(t, neg_2sig(1,:), 'r--');
    title('UKF Estimated State Trajectory for Given Data Log, \DeltaT=0.1sec');
    ylabel('\xi_g (m)');
    
    subplot(6,1,2);
    plot(t, state_est_errs(2,:), 'b'); hold on;
    plot(t, pos_2sig(2,:), 'r--'); hold on;
    plot(t, neg_2sig(2,:), 'r--');
    ylabel('\eta_g (m)');
    
    subplot(6,1,3);
    plot(t, wrapToPi(state_est_errs(3,:)), 'b'); hold on;
    plot(t, wrapToPi(state_est_errs(3,:))+sigs(3,:), 'r--'); hold on;
    plot(t, wrapToPi(state_est_errs(3,:))-sigs(3,:), 'r--');
    ylabel('\theta_g (rad)');
    
    subplot(6,1,4);
    plot(t, state_est_errs(4,:), 'b'); hold on;
    plot(t, pos_2sig(4,:), 'r--'); hold on;
    plot(t, neg_2sig(4,:), 'r--');
    ylabel('\xi_a (m)');
    
    subplot(6,1,5);
    plot(t, state_est_errs(5,:), 'b'); hold on;
    plot(t, pos_2sig(5,:), 'r--'); hold on;
    plot(t, neg_2sig(5,:), 'r--');
    ylabel('\eta_a (m)');
    
    subplot(6,1,6);
    plot(t, wrapToPi(state_est_errs(6,:)), 'b'); hold on;
    plot(t, wrapToPi(state_est_errs(6,:))+sigs(6,:), 'r--'); hold on;
    plot(t, wrapToPi(state_est_errs(6,:))-sigs(6,:), 'r--');
    xlabel('time (sec)'); ylabel('\theta_a (rad)');
end