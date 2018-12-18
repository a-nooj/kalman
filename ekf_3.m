rng(100)

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
        
P_plus_0 = 2.*[5 0 0 5 0 0; ...
			0 5 0 0 5 0; ...
			0 0 0.02 0 0 0; ...
			5 0 0 5 0 0; ...
			0 5 0 0 5 0; ...
			0 0 0 0 0 0.01];  % make smaller
        %which states uncertain? what expc vary tog etc

Qkf = 10.*[Qtrue(1,1) 0          0          0          0          0; ...
	   0          Qtrue(2,2) 0          0          0          0; ...
	   0          0          0.1*Qtrue(3,3) 0          0          0; ...
	   0          0          0          0.001*Qtrue(4,4) 0          0; ...
	   0          0          0          0          0.001*Qtrue(5,5) 0; ...
	   0          0          0          0          0          0.1*Qtrue(6,6)];
%P_plus_0=10.*eye(6,6);
%Qkf=Qtrue;
   
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
P_plus_0 = 1*[1 0 0 0 0 0; ...
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

x_hat_plus_prev = x_hat_plus_0;
P_plus_prev = P_plus_0;

x_hat_plus = [x_hat_plus_prev];
P_plus = [P_plus_prev];
sigs = [2.*sqrt(diag(P_plus_prev))];

for k = 1:length(tvec)-1
    t0 = tvec(k); t1 = tvec(k+1);

	% prediction
	q = zeros(num_states,1);  % simulating with w_k=0, function name below is misleading
	[x_hat_minus_k1, ~] = generate_noisy_ground_truth(u0, L, q, Rtrue, t0, t1, x_hat_plus_prev);
	[Ftilde_k, ~, ~, omegatilde_k] = get_discrete_mats(x_hat_plus_prev, u0, L, delta_T);
	P_minus_k1 = Ftilde_k*P_plus_prev*Ftilde_k' + omegatilde_k*Qkf*omegatilde_k';

	% measurement update
	y_hat_minus_k1 = get_measurement_value(x_hat_minus_k1);  % v_k+1 = 0
	y_hat_minus_k1(1) = wrapToPi(y_hat_minus_k1(1)); y_hat_minus_k1(3) = wrapToPi(y_hat_minus_k1(3));
    [~, ~, Htilde_k1, ~] = get_discrete_mats(x_hat_minus_k1, u0, L, delta_T);
	etilde_y_k1 = yreal(:,k+1) - y_hat_minus_k1;
    etilde_y_k1(1) = wrapToPi(etilde_y_k1(1)); etilde_y_k1(3) = wrapToPi(etilde_y_k1(3));
	S = Htilde_k1*P_minus_k1*Htilde_k1' + Rtrue;
	Ktilde_k1 = P_minus_k1*Htilde_k1'*inv(S);
	x_hat_plus_k1 = x_hat_minus_k1 + Ktilde_k1*etilde_y_k1;
	tmp = Ktilde_k1*Htilde_k1;
	P_plus_k1 = (eye(size(tmp)) - tmp)*P_minus_k1;
        
    % bookkeeping
    x_hat_plus_prev = x_hat_plus_k1;
    P_plus_prev = P_plus_k1;
    dx_hat_plus = [dx_hat_plus, dx_hat_plus_prev];
    P_plus = [P_plus, P_plus_prev];
    x_hat_plus = [x_hat_plus, x_hat_plus_prev];
    sigs = [sigs, 2.*sqrt(diag(P_plus_prev))];
end

pos_2sig = x_hat_plus + sigs; neg_2sig = x_hat_plus - sigs;

plot_state_errors(tvec, x_hat_plus, pos_2sig, neg_2sig,sigs);