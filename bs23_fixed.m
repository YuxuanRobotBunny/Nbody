function [t_hist, y_hist, stats] = bs23_fixed(rhs_fun, tspan, y0, h, params)
%BS23_FIXED Fixed-step Bogacki-Shampine RK23 integrator.
%
%   [t_hist, y_hist, stats] = BS23_FIXED(rhs_fun, [t0 tf], y0, h, params)
%   returns times and state history (columns correspond to times).

if numel(tspan) ~= 2
    error('tspan must be [t0 tf].');
end
if h <= 0
    error('Step size h must be positive.');
end

t0 = tspan(1);
tf = tspan(2);
if tf < t0
    error('This solver currently supports forward integration only.');
end

n_state = numel(y0);
n_est = max(2, ceil((tf - t0) / h) + 2);
t_hist = zeros(n_est, 1);
y_hist = zeros(n_state, n_est);
t_hist(1) = t0;
y_hist(:, 1) = y0;

n = 1;
nfev = 0;
t_eps = 1e-14 * max(1.0, abs(tf));

while t_hist(n) < tf - t_eps
    h_step = min(h, tf - t_hist(n));
    [y_next, ~, nfev_step] = bs23_step(rhs_fun, t_hist(n), y_hist(:, n), h_step, params);
    nfev = nfev + nfev_step;

    n = n + 1;
    if n > numel(t_hist)
        t_hist = [t_hist; zeros(n_est, 1)]; 
        y_hist = [y_hist, zeros(n_state, n_est)]; 
    end
    t_hist(n) = t_hist(n-1) + h_step;
    y_hist(:, n) = y_next;
end

t_hist = t_hist(1:n);
y_hist = y_hist(:, 1:n);

stats = struct();
stats.nfev = nfev;
stats.accepted_steps = n - 1;
stats.rejected_steps = 0;
stats.h_nominal = h;
end
