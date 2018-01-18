function dydt = thrd_ord_sys(t,y, ut, u, param)
% Third order system based on linear motor + msd sysem, Va input f output
% see IMECE paper - Bijo

q1 = param(1);
q2 = param(2);
q3 = param(3);
q4 = param(4);

% diff input
u_diff = [diff(u(1:2)); diff(u)];

% interpolate input
u_in = interp1(ut, u, t);
u_ind = interp1(ut, u_diff, t);

% Fddd = -q2/q1*fdd - q3/q1*fd - q4/q1*0 + q5/q1*Va
% equal to 
% Fddd = -q1*fdd - q2*fd + q3*Va
% this simplifies optimization -> fewer variables! 

dydt = [y(2);...
    y(3);...
    -q1*y(3) - q2*y(2)  + q3*u_ind + q4*u_in]; 


end