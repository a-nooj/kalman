function [] = plot_nees_nis(t, nees, nis, r1_nees, r2_nees, r1_nis, r2_nis)    
    figure;
    
    subplot(1,2,1);
    plot(t, repmat(r1_nees,[1,length(nees)]),'r','LineWidth',1.5); hold on;
    plot(t, repmat(r2_nees,[1,length(nees)]),'r','LineWidth',1.5); hold on;
    scatter(t, nees, 'b');
    title('NEES Estimation Results for EKF, 50 simulations'); xlabel('time step, k'); ylabel('NEES statistic, \bar{\epsilon}_x');
    legend('r_1 bound', 'r_2 bound', 'NEES @ time step k');
    
    subplot(1,2,2);
    plot(t, repmat(r1_nis,[1,length(nis)]),'r','LineWidth',1.5); hold on;
    plot(t, repmat(r2_nis,[1,length(nis)]),'r','LineWidth',1.5); hold on;
    scatter(t, nis, 'b');
    title('NIS Estimation Results for UKF, 50 simulations'); xlabel('time step, k'); ylabel('NIS statistic, \bar{\epsilon}_y');
    legend('r_1 bound', 'r_2 bound', 'NIS @ time step k');
end