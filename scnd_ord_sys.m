function dydt = scnd_ord_sys(t,y, ut, u, param)

% Standard 2nd Ord m-s-d system

q1 = param(1);
q2 = param(2);
q3 = param(3);

% interpolate input
u_in = interp1(ut, u, t);

dydt = [y(2); -q1*y(2) - q2*y(1) + u_in*q3]; % y' = y2

%q1 = b/m q2 = k/m q3 = 1/m
end