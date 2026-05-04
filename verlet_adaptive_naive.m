function [t_hist, R_hist, V_hist, stats] = verlet_adaptive_naive(accel_fun, tspan, R0, V0, h0, tol, params, opts)
%VERLET_ADAPTIVE_NAIVE Naive variable-step Verlet with step-doubling control.
%
%   This method is intentionally "non-geometric":
%   it adapts step size using local error estimates and therefore does not
%   preserve the symplectic structure in general.
%
%   [t_hist, R_hist, V_hist, stats] = VERLET_ADAPTIVE_NAIVE(accel_fun, tspan, R0, V0, h0, tol, params, opts)
%
%   opts fields (all optional):
%     rtol, atol         : scale for error ratio (defaults derived from tol)
%     safety             : step-size safety factor (default 0.9)
%     fac_min, fac_max   : min/max accepted growth factor (default 0.2, 2.0)
%     fac_max_reject     : max factor after rejection (default 0.8)
%     h_min, h_max       : step-size bounds
%     max_attempts       : hard cap on total attempts

if nargin < 8 || isempty(opts)
    opts = struct();
end
if numel(tspan) ~= 2
    error('tspan must be [t0 tf].');
end
if h0 <= 0
    error('Initial step size h0 must be positive.');
end
if tol <= 0
    error('tol must be positive.');
end
if ~isequal(size(R0), size(V0))
    error('R0 and V0 must have identical sizes.');
end

t0 = tspan(1);
tf = tspan(2);
if tf < t0
    error('This solver currently supports forward integration only.');
end

rtol = get_opt(opts, 'rtol', tol);
atol = get_opt(opts, 'atol', max(1e-12, 1e-3 * tol));
safety = get_opt(opts, 'safety', 0.9);
fac_min = get_opt(opts, 'fac_min', 0.2);
fac_max = get_opt(opts, 'fac_max', 2.0);
fac_max_reject = get_opt(opts, 'fac_max_reject', 0.8);
h_min = get_opt(opts, 'h_min', 1e-14 * max(1.0, abs(tf - t0)));
h_max = get_opt(opts, 'h_max', max(h0, tf - t0));
max_attempts = get_opt(opts, 'max_attempts', 2e6);

[N, d] = size(R0);
h = min(max(h0, h_min), h_max);
t_eps = 1e-14 * max(1.0, abs(tf));

n_est = max(2, ceil((tf - t0) / max(h, h_min)) + 2);
t_hist = zeros(n_est, 1);
R_hist = zeros(N, d, n_est);
V_hist = zeros(N, d, n_est);
h_accepted = zeros(n_est, 1);
err_accepted = zeros(n_est, 1);

t_hist(1) = t0;
R_hist(:, :, 1) = R0;
V_hist(:, :, 1) = V0;

attempt_cap = max(8, 2 * n_est);
t_attempt = zeros(attempt_cap, 1);
h_attempt = zeros(attempt_cap, 1);
err_attempt = zeros(attempt_cap, 1);
accepted_attempt = false(attempt_cap, 1);

A = accel_fun(R0, params);
if ~isequal(size(A), [N, d])
    error('accel_fun returned wrong size. Expected %d-by-%d.', N, d);
end
accel_evals = 1;

n = 1;
na = 0;

while t_hist(n) < tf - t_eps
    if na >= max_attempts
        error('Maximum attempted steps exceeded.');
    end
    if h < h_min
        error('Step size dropped below h_min.');
    end

    h_step = min(h, tf - t_hist(n));
    Rn = R_hist(:, :, n);
    Vn = V_hist(:, :, n);

    % One full step.
    [R_full, V_full, ~, eval1] = verlet_step_with_accel(accel_fun, Rn, Vn, A, h_step, params);
    accel_evals = accel_evals + eval1;

    % Two half steps.
    [R_half, V_half, A_half, eval2] = verlet_step_with_accel(accel_fun, Rn, Vn, A, 0.5 * h_step, params);
    [R_two, V_two, A_two, eval3] = verlet_step_with_accel(accel_fun, R_half, V_half, A_half, 0.5 * h_step, params);
    accel_evals = accel_evals + eval2 + eval3;

    err_ratio = position_error_ratio(R_two - R_full, Rn, R_two, atol, rtol);
    accept = err_ratio <= 1.0;

    na = na + 1;
    if na > attempt_cap
        grow = max(attempt_cap, 1024);
        t_attempt = [t_attempt; zeros(grow, 1)]; %
        h_attempt = [h_attempt; zeros(grow, 1)]; %
        err_attempt = [err_attempt; zeros(grow, 1)]; %
        accepted_attempt = [accepted_attempt; false(grow, 1)]; %
        attempt_cap = attempt_cap + grow;
    end
    t_attempt(na) = t_hist(n);
    h_attempt(na) = h_step;
    err_attempt(na) = err_ratio;
    accepted_attempt(na) = accept;

    if accept
        n = n + 1;
        if n > numel(t_hist)
            grow = max(numel(t_hist), 1024);
            t_hist = [t_hist; zeros(grow, 1)]; %
            R_hist = cat(3, R_hist, zeros(N, d, grow)); %
            V_hist = cat(3, V_hist, zeros(N, d, grow)); %
            h_accepted = [h_accepted; zeros(grow, 1)]; %
            err_accepted = [err_accepted; zeros(grow, 1)]; %
        end

        t_hist(n) = t_hist(n-1) + h_step;
        R_hist(:, :, n) = R_two;
        V_hist(:, :, n) = V_two;
        A = A_two;
        h_accepted(n-1) = h_step;
        err_accepted(n-1) = err_ratio;
    end

    h = next_h(h_step, err_ratio, accept, safety, fac_min, fac_max, fac_max_reject, h_min, h_max);
end

t_hist = t_hist(1:n);
R_hist = R_hist(:, :, 1:n);
V_hist = V_hist(:, :, 1:n);
if n > 1
    h_accepted = h_accepted(1:n-1);
    err_accepted = err_accepted(1:n-1);
else
    h_accepted = zeros(0, 1);
    err_accepted = zeros(0, 1);
end
t_attempt = t_attempt(1:na);
h_attempt = h_attempt(1:na);
err_attempt = err_attempt(1:na);
accepted_attempt = accepted_attempt(1:na);

stats = struct();
stats.accepted_steps = n - 1;
stats.rejected_steps = na - (n - 1);
stats.accel_evals = accel_evals;
stats.h_accepted = h_accepted;
stats.err_accepted = err_accepted;
stats.t_attempt = t_attempt;
stats.h_attempt = h_attempt;
stats.err_attempt = err_attempt;
stats.accepted_attempt = accepted_attempt;
stats.rtol = rtol;
stats.atol = atol;
end

function [R_next, V_next, A_next, evals] = verlet_step_with_accel(accel_fun, R, V, A, h, params)
R_next = R + h * V + 0.5 * h^2 * A;
A_next = accel_fun(R_next, params);
V_next = V + 0.5 * h * (A + A_next);
evals = 1;
end

function err_ratio = position_error_ratio(dR, R_old, R_new, atol, rtol)
scale = atol + rtol .* max(abs(R_old), abs(R_new));
scale(scale <= 0) = 1e-16;
err_ratio = max(abs(dR(:)) ./ scale(:));
end

function h_new = next_h(h, err_ratio, accept, safety, fac_min, fac_max, fac_max_reject, h_min, h_max)
if err_ratio <= 0
    err_ratio = 1e-16;
end
% Step-doubling estimator is O(h^3), so exponent = 1/3.
k = 1/3;
fac = safety * err_ratio^(-k);
if accept
    fac = min(fac_max, max(fac_min, fac));
else
    fac = min(fac_max_reject, max(fac_min, fac));
end
h_new = min(h_max, max(h_min, h * fac));
end

function v = get_opt(opts, name, default_value)
if isfield(opts, name)
    v = opts.(name);
else
    v = default_value;
end
end
