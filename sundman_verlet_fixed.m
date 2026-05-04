function [tau_hist, t_hist, R_hist, V_hist, stats] = sundman_verlet_fixed(accel_fun, g_fun, tspan, R0, V0, h_tau, params, opts)
%SUNDMAN_VERLET_FIXED Sundman-transformed leapfrog with fixed fictitious step.
%
%   This integrates
%       dR/dt = V, dV/dt = F(R)
%   using a time reparameterization
%       dt/dtau = g(R),  g(R) > 0
%   and a fixed step in tau.
%
%   The update is a leapfrog-like scheme in tau:
%       V_{n+1/2} = V_n + 0.5*dtau*g(R_n)*F(R_n)
%       R_{n+1}   = R_n + dtau*g(R_n)*V_{n+1/2}
%       t_{n+1}   = t_n + dtau*g(R_n)
%       V_{n+1}   = V_{n+1/2} + 0.5*dtau*g(R_{n+1})*F(R_{n+1})
%
%   Notes:
%   - This is Sundman-inspired and often behaves better than naive adaptive
%     stepping for eccentric motion.
%   - It is not a blanket guarantee of full symplecticity for arbitrary g.
%
%   Inputs:
%       accel_fun : A = accel_fun(R, params), size(A)=size(R)
%       g_fun     : g = g_fun(R, params), scalar > 0
%       tspan     : [t0 tf]
%       R0, V0    : N-by-d initial state
%       h_tau     : fixed nominal step in tau
%       params    : model parameters passed to accel_fun and g_fun
%       opts      : optional struct
%           opts.g_min, opts.g_max : clamp g(R)
%
%   Outputs:
%       tau_hist, t_hist : fictitious and physical time grids
%       R_hist, V_hist   : state history, size N-by-d-by-nt
%       stats            : counters and step diagnostics

if nargin < 8
    opts = struct();
end
if numel(tspan) ~= 2
    error('tspan must be [t0 tf].');
end
if h_tau <= 0
    error('h_tau must be positive.');
end
if ~isequal(size(R0), size(V0))
    error('R0 and V0 must have identical sizes.');
end

t0 = tspan(1);
tf = tspan(2);
if tf < t0
    error('This solver currently supports forward integration only.');
end

g_min = get_opt(opts, 'g_min', 1e-6);
g_max = get_opt(opts, 'g_max', inf);

[N, d] = size(R0);
n_est = max(2, ceil((tf - t0) / h_tau) + 2);

tau_hist = zeros(n_est, 1);
t_hist = zeros(n_est, 1);
R_hist = zeros(N, d, n_est);
V_hist = zeros(N, d, n_est);
dt_phys_hist = zeros(n_est, 1);

tau_hist(1) = 0.0;
t_hist(1) = t0;
R_hist(:, :, 1) = R0;
V_hist(:, :, 1) = V0;

A = accel_fun(R0, params);
if ~isequal(size(A), [N, d])
    error('accel_fun returned wrong size. Expected %d-by-%d.', N, d);
end
accel_evals = 1;

n = 1;
t_eps = 1e-14 * max(1.0, abs(tf));

while t_hist(n) < tf - t_eps
    Rn = R_hist(:, :, n);
    Vn = V_hist(:, :, n);

    g_n = clamp_g(g_fun(Rn, params), g_min, g_max);
    if ~(isscalar(g_n) && isfinite(g_n) && g_n > 0)
        error('g_fun must return a positive finite scalar.');
    end

    h_tau_step = min(h_tau, (tf - t_hist(n)) / g_n);
    if h_tau_step <= 0
        break;
    end

    V_half = Vn + 0.5 * h_tau_step * g_n * A;
    R_next = Rn + h_tau_step * g_n * V_half;
    t_next = t_hist(n) + h_tau_step * g_n;

    A_next = accel_fun(R_next, params);
    accel_evals = accel_evals + 1;
    g_next = clamp_g(g_fun(R_next, params), g_min, g_max);
    if ~(isscalar(g_next) && isfinite(g_next) && g_next > 0)
        error('g_fun must return a positive finite scalar.');
    end
    V_next = V_half + 0.5 * h_tau_step * g_next * A_next;

    n = n + 1;
    if n > numel(t_hist)
        grow = max(numel(t_hist), 1024);
        tau_hist = [tau_hist; zeros(grow, 1)]; %
        t_hist = [t_hist; zeros(grow, 1)]; %
        R_hist = cat(3, R_hist, zeros(N, d, grow)); %
        V_hist = cat(3, V_hist, zeros(N, d, grow)); %
        dt_phys_hist = [dt_phys_hist; zeros(grow, 1)]; %
    end

    tau_hist(n) = tau_hist(n-1) + h_tau_step;
    t_hist(n) = t_next;
    R_hist(:, :, n) = R_next;
    V_hist(:, :, n) = V_next;
    dt_phys_hist(n-1) = h_tau_step * g_n;

    A = A_next;
end

tau_hist = tau_hist(1:n);
t_hist = t_hist(1:n);
R_hist = R_hist(:, :, 1:n);
V_hist = V_hist(:, :, 1:n);
if n > 1
    dt_phys_hist = dt_phys_hist(1:n-1);
else
    dt_phys_hist = zeros(0, 1);
end

stats = struct();
stats.accepted_steps = n - 1;
stats.accel_evals = accel_evals;
stats.h_tau = h_tau;
stats.dt_phys_hist = dt_phys_hist;
if ~isempty(dt_phys_hist)
    stats.dt_phys_min = min(dt_phys_hist);
    stats.dt_phys_max = max(dt_phys_hist);
else
    stats.dt_phys_min = NaN;
    stats.dt_phys_max = NaN;
end
end

function g = clamp_g(g, g_min, g_max)
g = min(g_max, max(g_min, g));
end

function v = get_opt(opts, name, default_value)
if isfield(opts, name)
    v = opts.(name);
else
    v = default_value;
end
end
