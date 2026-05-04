function [t_hist, uhat_hist, stats] = sbdf2_kdv(grid, tspan, u0_hat, dt)
%SBDF2_KDV Semi-implicit BDF2 / AB2 scheme for KdV.

if numel(tspan) ~= 2
    error('tspan must be [t0 tf].');
end
if dt <= 0
    error('dt must be positive.');
end

t0 = tspan(1);
tf = tspan(2);
N = grid.N;
n_est = max(2, ceil((tf - t0) / dt) + 2);
t_hist = zeros(n_est, 1);
uhat_hist = zeros(N, n_est);
t_hist(1) = t0;
uhat_hist(:, 1) = u0_hat;

n = 1;
nonlinear_evals = 0;
t_eps = 1e-14 * max(1.0, abs(tf));

if tf > t0 + t_eps
    dt1 = min(dt, tf - t0);
    [t1, uhat1, stats1] = etdrk2_kdv(grid, [t0, t0 + dt1], u0_hat, dt1);
    n = 2;
    t_hist(2) = t1(end);
    uhat_hist(:, 2) = uhat1(:, end);
    nonlinear_evals = nonlinear_evals + stats1.nonlinear_evals;
end

while t_hist(n) < tf - t_eps
    dt_step = min(dt, tf - t_hist(n));
    Bn = kdv_nonlinear_hat(uhat_hist(:, n), grid);
    Bnm1 = kdv_nonlinear_hat(uhat_hist(:, n - 1), grid);
    nonlinear_evals = nonlinear_evals + 2;

    numer = 2 * uhat_hist(:, n) - 0.5 * uhat_hist(:, n - 1) + dt_step * (2 * Bn - Bnm1);
    denom = 1.5 - dt_step * grid.lin;
    uhat_next = numer ./ denom;
    uhat_next(grid.mean_idx) = uhat_hist(grid.mean_idx, 1);
    if grid.has_unmatched_mode
        uhat_next(grid.unmatched_idx) = 0.0;
    end

    n = n + 1;
    if n > numel(t_hist)
        grow = max(numel(t_hist), 1024);
        t_hist = [t_hist; zeros(grow, 1)]; %
        uhat_hist = [uhat_hist, zeros(N, grow)]; %
    end
    t_hist(n) = t_hist(n - 1) + dt_step;
    uhat_hist(:, n) = uhat_next;
end

t_hist = t_hist(1:n);
uhat_hist = uhat_hist(:, 1:n);
stats = struct('nonlinear_evals', nonlinear_evals, 'dt_nominal', dt, ...
    'method', 'SBDF2');
end
