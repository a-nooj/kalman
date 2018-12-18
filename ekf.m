%rng(100)

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
num_sims = 50;  %50

num_states = length(Qtrue);  % n
num_measurements = length(Rtrue);  % p

x_hat_plus_0 = x0;%[0;0;0;0;0;0];%x0;  % + [0; 1; 0; 0; 0; 0.1];
%{
P_plus_0 = 3.*[6 0 0 3 0 0; ...
			0 6 0 0 3 0; ...
			0 0 0.05 0 0 0; ...
			3 0 0 6 0 0; ...
			0 3 0 0 6 0; ...
			0 0 0 0 0 0.05];
%}
        
P_plus_0 = 2*[5 0 0 0 0 0; ...
			0 5 0 0 0 0; ...
			0 0 0.05 0 0 0; ...
			0 0 0 5 0 0; ...
			0 0 0 0 5 0; ...
			0 0 0 0 0 0.05];  % make smaller
        %which states uncertain? what expc vary tog etc

Qkf = 10.*[Qtrue(1,1) 0          0          0          0          0; ...
	   0          Qtrue(2,2) 0          0          0          0; ...
	   0          0          0.1*Qtrue(3,3) 0          0          0; ...
	   0          0          0          0.001*Qtrue(4,4) 0          0; ...
	   0          0          0          0          0.001*Qtrue(5,5) 0; ...
	   0          0          0          0          0          0.1*Qtrue(6,6)];
   
%{
P_plus_0 = 1.*[1 0 0 0 0 0; ...
			0 1 0 0 0 0; ...
			0 0 1 0 0 0; ...
			0 0 0 1 0 0; ...
			0 0 0 0 1 0; ...
			0 0 0 0 0 1];  % make smaller
        %which states uncertain? what expc vary tog etc

Qkf = 10.*[Qtrue(1,1) 0          0          0          0          0; ...
	   0          Qtrue(2,2) 0          0          0          0; ...
	   0          0          Qtrue(3,3) 0          0          0; ...
	   0          0          0          Qtrue(4,4) 0          0; ...
	   0          0          0          0          Qtrue(5,5) 0; ...
	   0          0          0          0          0          Qtrue(6,6)];
%}
%P_plus_0 = eye(6,6);
P_plus_0 = 2*[1 0 0 0 0 0; ...
			0 1 0 0 0 0; ...
			0 0 0.01 0 0 0; ...
			0 0 0 1 0 0; ...
			0 0 0 0 1 0; ...
			0 0 0 0 0 0.01];
Qkf=[Qtrue(1,1) 0          0          0          0          0; ...
	   0          Qtrue(2,2) 0          0          0          0; ...
	   0          0          Qtrue(3,3) 0          0          0; ...
	   0          0          0          Qtrue(4,4) 0          0; ...
	   0          0          0          0          Qtrue(5,5) 0; ...
	   0          0          0          0          0          Qtrue(6,6)];
Rkf = Rtrue;

nees = []; nis = [];

for i = 1:num_sims
	i
	x_hat_plus_prev = x_hat_plus_0;
	P_plus_prev = P_plus_0;

	x_hat_plus = [x_hat_plus_prev];
	P_plus = [P_plus_prev];
	sigs = [2.*sqrt(diag(P_plus_prev))];
	nees_trial = [nan]; nis_trial = [nan];

	xnoisy_prev = mvnrnd(x_hat_plus_prev, P_plus_prev)';
	noisy_state_traj = [xnoisy_prev];

	for k = 1:length(tvec)-1
        t0 = tvec(k); t1 = tvec(k+1);

        q = mvnrnd(zeros(1,num_states), Qtrue)';
		[xnoisy, ynoisy] = generate_noisy_ground_truth(u0, L, q, Rtrue, t0, t1, xnoisy_prev);
        
		% prediction
		q = zeros(num_states,1);  % simulating with w_k=0, function name below is misleading
		[x_hat_minus_k1, ~] = generate_noisy_ground_truth(u0, L, q, Rtrue, t0, t1, x_hat_plus_prev);
		[Ftilde_k, ~, ~, omegatilde_k] = get_discrete_mats(x_hat_plus_prev, u0, L, delta_T);
		P_minus_k1 = Ftilde_k*P_plus_prev*Ftilde_k' + omegatilde_k*Qkf*omegatilde_k';

		% measurement update
		y_hat_minus_k1 = get_measurement_value(x_hat_minus_k1);  % v_k+1 = 0
		y_hat_minus_k1(1) = wrapToPi(y_hat_minus_k1(1)); y_hat_minus_k1(3) = wrapToPi(y_hat_minus_k1(3));
        [~, ~, Htilde_k1, ~] = get_discrete_mats(x_hat_minus_k1, u0, L, delta_T);
		etilde_y_k1 = ynoisy - y_hat_minus_k1;
        etilde_y_k1(1) = wrapToPi(etilde_y_k1(1)); etilde_y_k1(3) = wrapToPi(etilde_y_k1(3));
		S = Htilde_k1*P_minus_k1*Htilde_k1' + Rtrue;
		Ktilde_k1 = P_minus_k1*Htilde_k1'*inv(S);
		x_hat_plus_k1 = x_hat_minus_k1 + Ktilde_k1*etilde_y_k1;
		tmp = Ktilde_k1*Htilde_k1;
		P_plus_k1 = (eye(size(tmp)) - tmp)*P_minus_k1;
            
        % bookkeeping
        x_hat_plus_prev = x_hat_plus_k1;
        P_plus_prev = P_plus_k1;
        xnoisy_prev = xnoisy;
		noisy_state_traj = [noisy_state_traj, xnoisy];
        dx_hat_plus = [dx_hat_plus, dx_hat_plus_prev];
        P_plus = [P_plus, P_plus_prev];
        x_hat_plus = [x_hat_plus, x_hat_plus_prev];
        sigs = [sigs, 2.*sqrt(diag(P_plus_prev))];

        % consistency tests
        x_k1 = xnoisy; e_x = x_k1 - x_hat_plus_k1;
        nees_stat = e_x'*inv(P_plus_k1)*e_x;
        nees_trial = [nees_trial, nees_stat];
        
        nis_stat = etilde_y_k1'*inv(S)*etilde_y_k1;
        nis_trial = [nis_trial, nis_stat];
	end

	state_est_errs = noisy_state_traj - x_hat_plus;
	pos_2sig = sigs; neg_2sig = -sigs;

	nees = [nees; nees_trial]; nis = [nis; nis_trial];
end

nees = mean(nees);
r1_nees = chi2inv(alpha_sig/2, num_sims*num_states)./num_sims;
r2_nees = chi2inv(1-(alpha_sig/2), num_sims*num_states)./num_sims;
percent_within_nees = 100*sum(nees>=r1_nees & nees<=r2_nees)/length(nees);

nis = mean(nis);
r1_nis = chi2inv(alpha_sig/2, num_sims*num_measurements)./num_sims;
r2_nis = chi2inv(1-(alpha_sig/2), num_sims*num_measurements)./num_sims;
percent_within_nis = 100*sum(nis>=r1_nis & nis<=r2_nis)/length(nis);

plot_nees_nis(tvec, nees, nis, r1_nees, r2_nees, r1_nis, r2_nis);
plot_state_errors(tvec, state_est_errs, pos_2sig, neg_2sig);