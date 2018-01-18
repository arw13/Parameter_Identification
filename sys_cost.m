function J_sysid = sys_cost(param, u, ut, y_exp)

ord = 2;

y0 = y_exp(1,1:ord);

if ord == 2
    [~,y_out] = ode45(@(t,y) scnd_ord_sys(t,y,ut,u, param), ut, y0);
elseif ord == 3
    [~,y_out] = ode45(@(t,y) thrd_ord_sys(t,y,ut,u, param), ut, y0);
end

m = length(y_out);

J_sysid = 1/m*sum(sum((y_exp(:,1:2) - y_out(:,1:2)).^2));