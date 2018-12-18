function dxdt = simulate_nonlinear_states(t, x, u, L, noisy, Q, noise_mat)
    dxdt = [u(1)*cos(x(3)); 
            u(1)*sin(x(3)); 
            u(1)*tan(u(2))/L; 
            u(3)*cos(x(6)); 
            u(3)*sin(x(6)); 
            u(4)];
    if noisy
        %tt = round(t);
        %if tt == 0
        %    tt=tt+1;
        %end
        %tt = floor(t/0.1)+1;
        %dxdt = dxdt + noise_mat(tt,:)';%mvnrnd(zeros(size(dxdt)), Q)';
        dxdt = dxdt + noise_mat;
    end
end