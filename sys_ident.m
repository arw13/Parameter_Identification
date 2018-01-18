%% SYS_Ident
clear

% Load training data and plot
filename = 'SYS_ID5.txt';
delimiterIn = ',';
headerlinesIn = 0;
A = importdata(filename,delimiterIn,headerlinesIn);

% Sys ID4 starts at sample 15
s = 4;
PWM = A(s:end,2);
Va = 12*PWM./256;
n = length(Va);
Ts = .01; %.1;
t = Ts*(1:n);

% F derivatives
F_meas = zeros(n,3);
F_meas(:,1) = A(s:end,1);
F_meas(:,2) = [mean(diff(F_meas(1:10,1))); diff(F_meas(:,1))]./Ts;
F_meas(:,3) = smooth([mean(diff(F_meas(1:2,2))); diff(F_meas(:,2))]./Ts, 3);

save('traindata.mat', 'F_meas', 'Va', 't', 'n')
%% optimization
% Call cost func from ode solver in fminsearch
% Constraints and initialization

init_params = [30 15 0.1];
% Best 2nd order fminsearch -> fit up to y(:,2), fval = 1609 init_params = [30 15 70];
% Best 2nd order patternsearch -> fit up to y(:,2), fval = 1607 init_params = [20 5 80];

% Best 3rd order fminsearch -> fit up to y(:,2), fval = 1600 init_params =
% [22 610 1400 500]; 1600    [40 700 1800 500]
% Best 3 term 3rd order fminsearch -> fit up to y(:,2), fval = 17910  init_params = [50 500 1500];

% Best updated 3rd order fmin -> fval = 2293 init_params =  [50 5 1.8e4
% 100]; opt_params = [58 5 13969 114]
% Best updated 3rd order pattern -> fval = 1681 init_params = [70 70 1e4 400]
% Best updated 3rd order PSO -> fval = 2215 init_params = [70 70 15000
% 400]; opt_params = [332 5 86835 1215]


A = []; b = []; Aeq = []; beq = [];
lb = [0.01 , .01, .01];
ub = [1e3, 1e3, 1e5];
nonlcon = [];

ut = t;
u = Va;
%
ObjFun = @(param) sys_cost(param, Va, t,F_meas);

% Nelder-Mead
opts = optimset('Display','iter','PlotFcns',@optimplotfval);
% [opt_params, fval] = fminsearch(@(param) sys_cost(param, Va, t,...
%     F_meas), init_params, opts);

% Constrained Interior Point Optimization
opts2 = optimoptions('fmincon','Display','iter','PlotFcns',@optimplotfval);
[opt_params, fval] = fmincon(@(param) sys_cost(param, Va, t,...
    F_meas), init_params, A, b, Aeq, beq, lb, ub, nonlcon, opts2);

% Adaptive Mesh Minimization
opts3 = psoptimset('Display','iter','PlotFcn',@psplotbestf);
% [opt_params, fval] = patternsearch(ObjFun,init_params, A,b,Aeq,beq,lb,ub,nonlcon,...
%     opts3);

% Particle Swarm Optimization
opts4 = optimoptions('particleswarm','Display','iter','PlotFcns',...
     @pswplotbestf, 'HybridFcn',@fmincon);
numvars = 3;
% [opt_params, fval] = particleswarm(ObjFun,numvars,lb,ub,opts4);


%% Compare to train data
ord = 2;

y0 = F_meas(1,1:ord) ; % [F_meas(1,1) -28 0];

if ord == 2
    [~,y_opt] = ode45(@(t,y) scnd_ord_sys(t,y,ut,u, opt_params), ut, y0);
elseif ord == 3
    [~,y_opt] = ode45(@(t,y) thrd_ord_sys(t,y,ut,u, opt_params), ut, y0);
end

% Plot

figure(2), clf
subplot(3,1,1)
grid on
% yyaxis right
plot(t, Va, 'k')
ylabel('Voltage Input, V')
title('Estimated System vs Training Data')

subplot(3,1,2:3)
hold on
grid on
% yyaxis left
plot(t, F_meas(:,1), 'b')
plot(t, F_meas(:,2), 'g')
plot(t, y_opt(:,1), 'b-.')
plot(t, y_opt(:,2), 'g-.')
ylabel('Force, N')

% subplot(3,1,3)
% hold on
% grid on
% plot(t, F_meas(:,3), 'c')
% plot(t, y_opt(:,3), 'c-.')
% xlabel('Time, s')
% legend('Experimental Data', 'Estimated Model Response')

%% Compare to test data
filename = 'SYS_ID4.txt';
delimiterIn = ',';
headerlinesIn = 0;
A_test = importdata(filename,delimiterIn,headerlinesIn);

% Sys ID3 starts at sample 24
s = 4;
en = 39;
PWM_test = A_test(s:end,2);
Va_test = 12*PWM_test./256;
n_test = length(Va_test);
t_test = Ts*(1:n_test);
ut_test = t_test;
u_test = Va_test;

F_meas_test = zeros(n_test,3);
F_meas_test(:,1) = A_test(s:end,1);
F_meas_test(:,2) = [mean(diff(F_meas_test(1:10,1))); diff(F_meas_test(:,1))]./Ts;
F_meas_test(:,3) = [mean(diff(F_meas_test(1:10,2))); diff(F_meas_test(:,2))]./Ts;

% % test with params
y0_test = F_meas_test(1,1:ord); %[F_meas_test(1,1) -25 0];

if ord == 2
    [~,y_opt_test] = ode45(@(t,y) scnd_ord_sys(t,y,ut_test,u_test,...
        opt_params), ut_test, y0_test);
elseif ord ==3
    [~,y_opt_test] = ode45(@(t,y) thrd_ord_sys(t,y,ut_test,u_test,...
        opt_params), ut_test, y0_test);
end

figure(3), clf
subplot(3,1,1)
grid on
% yyaxis right
plot(t_test, Va_test, 'k')
ylabel('Voltage Input, V')
title('Estimated System vs Test Data')

subplot(3,1,2:3)
hold on
grid on
% yyaxis left
plot(t_test, F_meas_test(:,1), 'b') %force
plot(t_test, F_meas_test(:,2), 'g') %force deriv
% plot(t_test, F_meas_test(:,3), 'r-.')
plot(t_test, y_opt_test(:,1), 'b-.')
plot(t_test, y_opt_test(:,2), 'g-.')
% plot(t_test, y_opt_test(:,3), 'r')

ylabel('Force, N')
xlabel('Time, s')
% legend('Experimental Data', 'Estimated Model Response', 'Location', 'Northwest')

