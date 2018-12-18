function [] = get_plots(t, states, measurements)
    f1 = figure;
    f2 = figure;
    
    figure(f1);
    plot_states(t, states);
    
    figure(f2);
    a = load('cooplocalization_finalproj_KFdata.mat');
    plot_measurements(t, measurements, a.measLabels);
end

function [] = plot_states(t, states)
    subplot(6,1,1);
    plot(t,states(1,:));
    title('Noisy Simulated Ground Truth States vs. Time, \DeltaT=0.1 sec');
    ylabel('\xi_g (m)');

    subplot(6,1,2);
    plot(t,states(2,:));
    ylabel('\eta_g (m)');

    subplot(6,1,3);
    plot(t,wrapToPi(states(3,:)));
    ylabel('\theta_g (rad)');

    subplot(6,1,4);
    plot(t,states(4,:));
    ylabel('\xi_a (m)');

    subplot(6,1,5);
    plot(t,states(5,:));
    ylabel('\eta_a (m)');

    subplot(6,1,6);
    plot(t,wrapToPi(states(6,:)));
    xlabel('time (sec)'); ylabel('\theta_a (rad)');
end

function [] = plot_measurements(t, measurements, labels)
    subplot(5,1,1);
    plot(t,measurements(1,:));
    title('Noisy Simulated Data vs. Time, \DeltaT=0.1 sec');
    ylabel(labels(1));

    subplot(5,1,2);
    plot(t,measurements(2,:));
    ylabel(labels(2));

    subplot(5,1,3);
    plot(t,measurements(3,:));
    ylabel(labels(3));

    subplot(5,1,4);
    plot(t,measurements(4,:));
    ylabel(labels(4));

    subplot(5,1,5);
    plot(t,measurements(5,:));
    xlabel('time (sec)'); ylabel(labels(5));
end