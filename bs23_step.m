function [y3, err_vec, nfev] = bs23_step(rhs_fun, t, y, h, params)
%BS23_STEP One Bogacki-Shampine RK23 step.
%
%   y3      : third-order solution at t+h
%   err_vec : y3 - y2, where y2 is the embedded second-order estimate
%   nfev    : number of RHS evaluations (always 4)

k1 = rhs_fun(t, y, params);
k2 = rhs_fun(t + 0.5 * h, y + h * (0.5 * k1), params);
k3 = rhs_fun(t + 0.75 * h, y + h * (0.75 * k2), params);

y3 = y + h * ((2/9) * k1 + (1/3) * k2 + (4/9) * k3);
k4 = rhs_fun(t + h, y3, params);

y2 = y + h * ((7/24) * k1 + 0.25 * k2 + (1/3) * k3 + 0.125 * k4);
err_vec = y3 - y2;
nfev = 4;
end
