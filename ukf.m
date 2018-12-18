%rng(100)

eps=0.001;

delta_T = 0.1;
L = 0.5;
x0 = [10; 0; pi/2; -60; 0; -pi/2];
u0 = [2; -pi/18; 12; pi/25];

given = load('cooplocalization_finalproj_KFdata.mat');
Qtrue = given.Qtrue;
Rtrue = given.Rtrue;
yreal = given.ydata;
tvec = given.tvec;
measLabels = given.measLabels;

alpha_sig = 0.05;
num_sims = 1;  %50

num_states = length(Qtrue);  % n
num_measurements = length(Rtrue);  % p

x_hat_plus_0 = x0;%[0;0;0;0;0;0];%x0;  % + [0; 1; 0; 0; 0; 0.1];
P_plus_0 = [50 0 0 0 0 0; ...
			0 50 0 0 0 0; ...
			0 0 0.5 0 0 0; ...
			0 0 0 50 0 0; ...
			0 0 0 0 50 0; ...
			0 0 0 0 0 0.5];

Qkf = 0.03.*[Qtrue(1,1) 0          0          0          0          0; ...
	   0          Qtrue(2,2) 0          0          0          0; ...
	   0          0          Qtrue(3,3) 0          0          0; ...
	   0          0          0          Qtrue(4,4) 0          0; ...
	   0          0          0          0          Qtrue(5,5) 0; ...
	   0          0          0          0          0          Qtrue(6,6)];

 P_plus_0 = 5.*[2 0 0 0 0 0; ...
			0 2 0 0 0 0; ...
			0 0 0.1 0 0 0; ...
			0 0 0 2 0 0; ...
			0 0 0 0 2 0; ...
			0 0 0 0 0 0.1];  % make smaller
        %which states uncertain? what expc vary tog etc

Qkf = 0.008.*[Qtrue(1,1) 0          0          0          0          0; ...
	   0          Qtrue(2,2) 0          0          0          0; ...
	   0          0          Qtrue(3,3) 0          0          0; ...
	   0          0          0          0.0001*Qtrue(4,4) 0          0; ...
	   0          0          0          0          0.0001*Qtrue(5,5) 0; ...
	   0          0          0          0          0          Qtrue(6,6)];

P_plus_0 = 2*[1 0 0 0 0 0; ...
			0 1 0 0 0 0; ...
			0 0 0.01 0 0 0; ...
			0 0 0 1 0 0; ...
			0 0 0 0 1 0; ...
			0 0 0 0 0 0.01];
Qkf = 1.*[Qtrue(1,1) 0          0          0          0          0; ...
	   0          Qtrue(2,2) 0          0          0          0; ...
	   0          0          Qtrue(3,3) 0          0          0; ...
	   0          0          0          Qtrue(4,4) 0          0; ...
	   0          0          0          0          Qtrue(5,5) 0; ...
	   0          0          0          0          0          Qtrue(6,6)];
   
Rkf = Rtrue;

nees = []; nis = [];

n = num_states;
alpha = 1; beta = 2; kappa = 0;
lambda = alpha^2 * (n + kappa) - n;

for nsim = 1:num_sims
	nsim
	x_hat_plus_prev = x_hat_plus_0;
	P_plus_prev = P_plus_0;

	x_hat_plus = [x_hat_plus_prev];
	P_plus = [P_plus_prev];
	sigs = [2.*sqrt(diag(P_plus_prev))];
	nees_trial = [nan]; nis_trial = [nan];

	%xnoisy_prev = mvnrnd(x_hat_plus_prev, P_plus_prev)';
	%noisy_state_traj = [xnoisy_prev];

	for k = 1:length(tvec)-1
        t0 = tvec(k); t1 = tvec(k+1);

        %q = mvnrnd(zeros(1,num_states), Qtrue)';
		%[xnoisy, ynoisy] = generate_noisy_ground_truth(u0, L, q, Rtrue, t0, t1, xnoisy_prev);

		% dynamics prediction
        S_k = chol(P_plus_prev);
        chi_0_k = x_hat_plus_prev;
		chi_i_k = [chi_0_k];
		for j=1:n
			chi_i_k_val = x_hat_plus_prev + sqrt(n + lambda)*S_k(j,:)';
			chi_i_k = [chi_i_k, chi_i_k_val];
		end
		for j=1:n
			chi_i_k_val = x_hat_plus_prev - sqrt(n + lambda)*S_k(j,:)';
			chi_i_k = [chi_i_k, chi_i_k_val];
		end

		chibar_i_k1 = [];
		for i = 1:size(chi_i_k, 2)
			q = zeros(num_states,1);  % simulating with w_k=0, function name below is misleading
			[chibar_i_k1_val, ~] = generate_noisy_ground_truth(u0, L, q, Rtrue, t0, t1, chi_i_k(:,i));
			chibar_i_k1 = [chibar_i_k1, chibar_i_k1_val];
		end

		w_m_0 = lambda/(n + lambda);
		w_c_0 = (lambda/(n + lambda)) + 1 - alpha^2 + beta;
		w_m_i = 1/(2*(n + lambda));
		w_c_i = w_m_i;
		x_hat_minus_k1 = zeros(num_states, 1);
		for i = 1:size(chibar_i_k1,2)
			if i == 1
				multiplier_x = w_m_0;
			else
				multiplier_x = w_m_i;
			end
			x_hat_minus_k1 = x_hat_minus_k1 + multiplier_x*chibar_i_k1(:,i);
        end
		P_minus_k1 = zeros(num_states, num_states);
        for i = 1:size(chibar_i_k1,2)
            if i == 1
				multiplier_P = w_c_0;
			else
				multiplier_P = w_c_i;
			end
            P_minus_k1 = P_minus_k1 + multiplier_P*(chibar_i_k1(:,i) - x_hat_minus_k1)*(chibar_i_k1(:,i) - x_hat_minus_k1)';  % Q_k outside?
        end
        P_minus_k1 = P_minus_k1 + Qkf;

		% measurement update
        Sbar_k1=chol(P_minus_k1);
		chi_0_k1 = x_hat_minus_k1;
		chi_i_k1 = [chi_0_k1];
		for j=1:n
			chi_i_k_val = x_hat_minus_k1 + sqrt(n + lambda)*Sbar_k1(j,:)';
			chi_i_k1 = [chi_i_k1, chi_i_k_val];
		end
		for j=1:n
			chi_i_k_val = x_hat_minus_k1 - sqrt(n + lambda)*Sbar_k1(j,:)';
			chi_i_k1 = [chi_i_k1, chi_i_k_val];
		end

		gamma_i_k1 = [];
		for i = 1:size(chi_i_k, 2)
			gamma_i_k1_val = get_measurement_value(chi_i_k1(:,i));
			gamma_i_k1 = [gamma_i_k1, gamma_i_k1_val];
		end

		y_hat_minus_k1 = zeros(num_measurements, 1);
		for i = 1:size(gamma_i_k1,2)
			if i == 1
				multiplier_y = w_m_0;
			else
				multiplier_y = w_m_i;
			end
			y_hat_minus_k1 = y_hat_minus_k1 + multiplier_y*gamma_i_k1(:,i);
        end
        %y_hat_minus_k1(1) = wrapToPi(y_hat_minus_k1(1)); y_hat_minus_k1(3) = wrapToPi(y_hat_minus_k1(3));
		P_yy_k1 = zeros(num_measurements, num_measurements);
        for i = 1:size(gamma_i_k1,2)
			if i == 1
				multiplier_P = w_c_0;
			else
				multiplier_P = w_c_i;
			end
			P_yy_k1 = P_yy_k1 + multiplier_P*(gamma_i_k1(:,i) - y_hat_minus_k1)*(gamma_i_k1(:,i) - y_hat_minus_k1)';  % R_k1 outside?
        end
        P_yy_k1=P_yy_k1 + Rkf;
        
		P_xy_k1 = zeros(num_states, num_measurements);
		for i = 1:size(chibar_i_k1,2)
			if i == 1
				multiplier_P = w_c_0;
			else
				multiplier_P = w_c_i;
			end
			P_xy_k1 = P_xy_k1 + multiplier_P*(chi_i_k1(:,i) - x_hat_minus_k1)*(gamma_i_k1(:,i) - y_hat_minus_k1)';  % slides say gamma_i_k...
		end

		K_k1 = P_xy_k1*inv(P_yy_k1);

        y_hat_minus_k1(1) = wrapToPi(y_hat_minus_k1(1)); y_hat_minus_k1(3) = wrapToPi(y_hat_minus_k1(3));
		e_y = yreal(:,k+1) - y_hat_minus_k1;
		e_y(1) = wrapToPi(e_y(1)); e_y(3) = wrapToPi(e_y(3));
		x_hat_plus_k1 = x_hat_minus_k1 + K_k1*e_y;
		P_plus_k1 = P_minus_k1 - K_k1*P_yy_k1*K_k1';  % alt formulation...
        
		% bookkeeping
        x_hat_plus_prev = x_hat_plus_k1;
        P_plus_prev = P_plus_k1;
        %xnoisy_prev = xnoisy;
		%noisy_state_traj = [noisy_state_traj, xnoisy];
        P_plus = [P_plus, P_plus_prev];
        x_hat_plus = [x_hat_plus, x_hat_plus_prev];
        sigs = [sigs, 2.*sqrt(diag(P_plus_prev))];

        % consistency tests
        %x_k1 = xnoisy; e_x = x_k1 - x_hat_plus_k1;
        %nees_stat = e_x'*inv(P_plus_k1)*e_x;
        %nees_trial = [nees_trial, nees_stat];
        
        nis_stat = e_y'*inv(P_yy_k1)*e_y;
        nis_trial = [nis_trial, nis_stat];
    end

	%state_est_errs = noisy_state_traj - x_hat_plus;
	%pos_2sig = sigs; neg_2sig = -sigs;

	nees = [nees; nees_trial]; nis = [nis; nis_trial];
end

%nees = mean(nees);
%r1_nees = chi2inv(alpha_sig/2, num_sims*num_states)./num_sims;
%r2_nees = chi2inv(1-(alpha_sig/2), num_sims*num_states)./num_sims;
%percent_within_nees = 100*sum(nees>=r1_nees & nees<=r2_nees)/length(nees);

%nis = mean(nis);
r1_nis = chi2inv(alpha_sig/2, num_sims*num_measurements)./num_sims;
r2_nis = chi2inv(1-(alpha_sig/2), num_sims*num_measurements)./num_sims;
percent_within_nis = 100*sum(nis_trial>=r1_nis & nis_trial<=r2_nis)/length(nis_trial);

plot_nees_nis(tvec*10, zeros(size(tvec)), nis_trial, 2, 4, r1_nis, r2_nis);
%plot_state_errors(tvec, state_est_errs, pos_2sig, neg_2sig);