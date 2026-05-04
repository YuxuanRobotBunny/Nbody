function [t_hist, R_hist, V_hist, stats] = verlet_fixed(accel_fun, tspan, R0, V0, h, params)
%VERLET_FIXED Fixed-step velocity-Verlet integrator (2nd-order symplectic).
%
%   [t_hist, R_hist, V_hist, stats] = VERLET_FIXED(accel_fun, [t0 tf], R0, V0, h, params)
%   integrates second-order ODEs in the form:
%       dR/dt = V, dV/dt = accel_fun(R, params)
%
%   Inputs:
%       accel_fun : function handle, A = accel_fun(R, params), size(A)=size(R)
%       tspan     : [t0 tf]
%       R0, V0    : N-by-d initial position and velocity
%       h         : fixed nominal step size
%       params    : parameter struct passed to accel_fun
%
%   Outputs:
%       t_hist    : time samples, length nt
%       R_hist    : N-by-d-by-nt positions
%       V_hist    : N-by-d-by-nt velocities
%       stats     : struct with step and evaluation counters

if numel(tspan) ~= 2
    error('tspan must be [t0 tf].');
end
if h <= 0
    error('Step size h must be positive.');
end
if ~isequal(size(R0), size(V0))
    error('R0 and V0 must have identical sizes.');
end

t0 = tspan(1);
tf = tspan(2);
if tf < t0
    error('This solver currently supports forward integration only.');
end

[N, d] = size(R0);
n_est = max(2, ceil((tf - t0) / h) + 2);

t_hist = zeros(n_est, 1);
R_hist = zeros(N, d, n_est);
V_hist = zeros(N, d, n_est);

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
    h_step = min(h, tf - t_hist(n));
    Rn = R_hist(:, :, n);
    Vn = V_hist(:, :, n);

    R_next = Rn + h_step * Vn + 0.5 * h_step^2 * A;
    A_next = accel_fun(R_next, params);
    accel_evals = accel_evals + 1;
    V_next = Vn + 0.5 * h_step * (A + A_next);

    n = n + 1;
    if n > numel(t_hist)
        grow = max(numel(t_hist), 1024);
        t_hist = [t_hist; zeros(grow, 1)]; %
        R_hist = cat(3, R_hist, zeros(N, d, grow)); %
        V_hist = cat(3, V_hist, zeros(N, d, grow)); %
    end

    t_hist(n) = t_hist(n-1) + h_step;
    R_hist(:, :, n) = R_next;
    V_hist(:, :, n) = V_next;
    A = A_next;
end

t_hist = t_hist(1:n);
R_hist = R_hist(:, :, 1:n);
V_hist = V_hist(:, :, 1:n);

stats = struct();
stats.accepted_steps = n - 1;
stats.rejected_steps = 0;
stats.h_nominal = h;
stats.accel_evals = accel_evals;
end
