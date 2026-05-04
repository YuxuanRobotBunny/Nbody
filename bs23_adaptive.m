function [t_hist, y_hist, stats] = bs23_adaptive(rhs_fun, tspan, y0, h0, tol, params, opts)
%BS23_ADAPTIVE Adaptive Bogacki-Shampine RK23 integrator.
%
%   [t_hist, y_hist, stats] = BS23_ADAPTIVE(rhs_fun, [t0 tf], y0, h0, tol, params, opts)
%   adapts h to keep local error within tolerance.
%
%   opts.err_mode:
%       'position' -> only first N*dim components are controlled
%       'full' (default) -> control all state components
%   opts.rtol / opts.atol:
%       weighted error scale uses atol + rtol*max(|y_n|, |y_{n+1}|) on controlled
%       components. If omitted, defaults are derived from tol:
%           rtol = tol
%           atol = max(1e-12, 1e-3*tol)

if nargin < 6
    params = [];
end
if nargin < 7 || isempty(opts)
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

t0 = tspan(1);
tf = tspan(2);
if tf < t0
    error('This solver currently supports forward integration only.');
end

err_mode = get_opt(opts, 'err_mode', 'full');
safety = get_opt(opts, 'safety', 0.9);
fac_min = get_opt(opts, 'fac_min', 0.2);
fac_max = get_opt(opts, 'fac_max', 5.0);
h_min = get_opt(opts, 'h_min', 1e-14 * max(1.0, abs(tf - t0)));
h_max = get_opt(opts, 'h_max', max(h0, tf - t0));
max_attempts = get_opt(opts, 'max_attempts', 2e6);
rtol = get_opt(opts, 'rtol', tol);
atol = get_opt(opts, 'atol', max(1e-12, 1e-3 * tol));
controller = get_opt(opts, 'controller', 'pi');
fac_max_reject = get_opt(opts, 'fac_max_reject', 0.9);

n_state = numel(y0);
t_eps = 1e-14 * max(1.0, abs(tf));
h = min(max(h0, h_min), h_max);

n_est = max(2, ceil((tf - t0) / max(h, h_min)) + 2);
t_hist = zeros(n_est, 1);
y_hist = zeros(n_state, n_est);
h_accepted = zeros(n_est, 1);
err_accepted = zeros(n_est, 1);

t_hist(1) = t0;
y_hist(:, 1) = y0;

attempt_cap = max(8, 2 * n_est);
t_attempt = zeros(attempt_cap, 1);
h_attempt = zeros(attempt_cap, 1);
err_attempt = zeros(attempt_cap, 1);
accepted_attempt = false(attempt_cap, 1);

n = 1;
na = 0;
nfev = 0;
err_prev_accept = NaN;

while t_hist(n) < tf - t_eps
    if na >= max_attempts
        error('Maximum number of attempted steps exceeded.');
    end
    if h < h_min
        error('Step size dropped below h_min. Likely near-collision or overly strict tol.');
    end

    h_step = min(h, tf - t_hist(n));
    [y3, err_vec, nfev_step] = bs23_step(rhs_fun, t_hist(n), y_hist(:, n), h_step, params);
    nfev = nfev + nfev_step;

    err_ratio = compute_error_ratio(err_vec, y_hist(:, n), y3, err_mode, params, atol, rtol);
    accept = err_ratio <= 1.0;

    na = na + 1;
    if na > attempt_cap
        grow = max(attempt_cap, 1024);
        t_attempt = [t_attempt; zeros(grow, 1)]; 
        h_attempt = [h_attempt; zeros(grow, 1)]; 
        err_attempt = [err_attempt; zeros(grow, 1)]; 
        accepted_attempt = [accepted_attempt; false(grow, 1)]; 
        attempt_cap = attempt_cap + grow;
    end
    t_attempt(na) = t_hist(n);
    h_attempt(na) = h_step;
    err_attempt(na) = err_ratio;
    accepted_attempt(na) = accept;

    err_last_accept = err_prev_accept;
    if accept
        n = n + 1;
        if n > numel(t_hist)
            grow = max(numel(t_hist), 1024);
            t_hist = [t_hist; zeros(grow, 1)]; 
            y_hist = [y_hist, zeros(n_state, grow)]; 
            h_accepted = [h_accepted; zeros(grow, 1)]; 
            err_accepted = [err_accepted; zeros(grow, 1)]; %
        end

        t_hist(n) = t_hist(n-1) + h_step;
        y_hist(:, n) = y3;
        h_accepted(n-1) = h_step;
        err_accepted(n-1) = err_ratio;
        err_prev_accept = max(err_ratio, 1e-16);
    end

    h = next_step_size(h_step, err_ratio, err_last_accept, accept, controller, ...
        safety, fac_min, fac_max, fac_max_reject, h_min, h_max);
end

t_hist = t_hist(1:n);
y_hist = y_hist(:, 1:n);
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
stats.nfev = nfev;
stats.accepted_steps = n - 1;
stats.rejected_steps = na - (n - 1);
stats.h_accepted = h_accepted;
stats.err_accepted = err_accepted;
stats.t_attempt = t_attempt;
stats.h_attempt = h_attempt;
stats.err_attempt = err_attempt;
stats.accepted_attempt = accepted_attempt;
stats.rtol = rtol;
stats.atol = atol;
stats.err_mode = err_mode;
end

function v = get_opt(opts, name, default_value)
if isfield(opts, name)
    v = opts.(name);
else
    v = default_value;
end
end

function h_new = next_step_size(h, err_ratio, err_prev_accept, accept, controller, ...
    safety, fac_min, fac_max, fac_max_reject, h_min, h_max)
if err_ratio <= 0
    err_ratio = 1e-16;
end

% Embedded error is order p=2, so exponent is 1/(p+1)=1/3.
k = 1/3;
alpha = 0.7 * k;
beta = 0.4 * k;

if accept
    if strcmpi(controller, 'pi') && isfinite(err_prev_accept) && err_prev_accept > 0
        fac = safety * err_ratio^(-alpha) * err_prev_accept^(beta);
    else
        fac = safety * err_ratio^(-k);
    end
    fac = min(fac_max, max(fac_min, fac));
else
    fac = safety * err_ratio^(-k);
    fac = min(fac_max_reject, max(fac_min, fac));
end
h_new = min(h_max, max(h_min, h * fac));
end

function err_ratio = compute_error_ratio(err_vec, y_old, y_new, err_mode, params, atol, rtol)
[idx, n_state] = controlled_indices(err_mode, params, numel(err_vec));
err_ctrl = abs(err_vec(idx));

atol_ctrl = expand_tolerance(atol, n_state, idx);
rtol_ctrl = expand_tolerance(rtol, n_state, idx);
scale = atol_ctrl + rtol_ctrl .* max(abs(y_old(idx)), abs(y_new(idx)));
scale(scale <= 0) = 1e-16;

err_ratio = max(err_ctrl ./ scale);
end

function [idx, n_state] = controlled_indices(err_mode, params, n_state)
switch lower(err_mode)
    case 'position'
        if ~isstruct(params) || ~isfield(params, 'N') || ~isfield(params, 'dim')
            error('err_mode ''position'' requires params.N and params.dim.');
        end
        npos = params.N * params.dim;
        idx = (1:npos).';
    case 'full'
        idx = (1:n_state).';
    otherwise
        error('Unknown err_mode: %s', err_mode);
end
end

function t_ctrl = expand_tolerance(t, n_state, idx)
if isscalar(t)
    t_ctrl = repmat(t, numel(idx), 1);
    return;
end
t = t(:);
if numel(t) == n_state
    t_ctrl = t(idx);
    return;
end
if numel(t) == numel(idx)
    t_ctrl = t;
    return;
end
error('Tolerance size mismatch. Use scalar, full-state vector, or controlled-subset vector.');
end
