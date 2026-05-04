function [tau_hist, t_hist, R_hist, V_hist, stats] = canonical_sundman_midpoint( ...
    mu, g_fun, grad_g_fun, tspan, R0, V0, h_tau, opts)
%CANONICAL_SUNDMAN_MIDPOINT
% Symplectic extended-phase Sundman integrator via implicit midpoint.
%
%   The physical Hamiltonian is:
%       H(R, V) = 0.5*|V|^2 - mu/|R|
%
%   Define an extended Hamiltonian:
%       K = g(R) * ( H(R,V) + p_t )
%
%   with fictitious time tau and canonical variables (R,V,t,p_t).
%   We integrate Hamilton's equations for K using implicit midpoint, which is
%   symplectic for Hamiltonian systems.
%
%   Inputs:
%     mu        : gravitational parameter (G*m0), scalar
%     g_fun     : g = g_fun(R), scalar > 0
%     grad_g_fun: grad g wrt R, column or row vector length d
%     tspan     : [t0 tf] in physical time
%     R0, V0    : initial position and velocity (1-by-d or d-by-1)
%     h_tau     : fixed nominal step in fictitious time
%     opts      : optional struct
%       opts.newton_tol
%       opts.newton_maxit
%       opts.g_min, opts.g_max (clamp g)
%
%   Outputs:
%     tau_hist, t_hist : fictitious and physical time histories
%     R_hist, V_hist   : 1-by-d-by-nt state history
%     stats            : diagnostics including constraint C=H+p_t

if nargin < 8
    opts = struct();
end
if numel(tspan) ~= 2
    error('tspan must be [t0 tf].');
end
if h_tau <= 0
    error('h_tau must be positive.');
end

R0 = R0(:);
V0 = V0(:);
if numel(R0) ~= numel(V0)
    error('R0 and V0 must have identical lengths.');
end

t0 = tspan(1);
tf = tspan(2);
if tf < t0
    error('This solver currently supports forward integration only.');
end

d = numel(R0);
nvar = 2 * d + 2;  % [R; V; t; p_t]

newton_tol = get_opt(opts, 'newton_tol', 1e-12);
newton_maxit = get_opt(opts, 'newton_maxit', 20);
g_min = get_opt(opts, 'g_min', 1e-8);
g_max = get_opt(opts, 'g_max', inf);

H0 = 0.5 * sum(V0.^2) - mu / norm(R0, 2);
pt0 = -H0;  % start on constraint manifold H + p_t = 0

x = [R0; V0; t0; pt0];

n_est = max(2, ceil((tf - t0) / h_tau) + 2);
tau_hist = zeros(n_est, 1);
t_hist = zeros(n_est, 1);
R_hist = zeros(1, d, n_est);
V_hist = zeros(1, d, n_est);
pt_hist = zeros(n_est, 1);
C_hist = zeros(n_est, 1);
dt_phys_hist = zeros(n_est, 1);
newton_it_hist = zeros(n_est, 1);
newton_res_hist = zeros(n_est, 1);

tau_hist(1) = 0.0;
t_hist(1) = t0;
R_hist(1, :, 1) = R0.';
V_hist(1, :, 1) = V0.';
pt_hist(1) = pt0;
C_hist(1) = H0 + pt0;

n = 1;
rhs_evals = 0;
t_eps = 1e-14 * max(1.0, abs(tf));

while x(2*d + 1) < tf - t_eps
    Rn = x(1:d);
    g_n = clamp_scalar(g_fun(Rn), g_min, g_max);
    if ~(isscalar(g_n) && isfinite(g_n) && g_n > 0)
        error('g_fun must return positive finite scalar.');
    end

    % Choose final step so physical time does not overshoot too much.
    h_step = min(h_tau, (tf - x(2*d + 1)) / g_n);
    if h_step <= 0
        break;
    end

    [x_next, nit, fres, ok, evals_step] = midpoint_step_newton( ...
        x, h_step, mu, g_fun, grad_g_fun, g_min, g_max, newton_tol, newton_maxit);
    rhs_evals = rhs_evals + evals_step;
    if ~ok
        error('Newton failed at step %d (residual %.3e).', n, fres);
    end

    n = n + 1;
    if n > numel(t_hist)
        grow = max(numel(t_hist), 1024);
        tau_hist = [tau_hist; zeros(grow, 1)]; 
        t_hist = [t_hist; zeros(grow, 1)]; 
        R_hist = cat(3, R_hist, zeros(1, d, grow)); 
        V_hist = cat(3, V_hist, zeros(1, d, grow)); 
        pt_hist = [pt_hist; zeros(grow, 1)]; 
        C_hist = [C_hist; zeros(grow, 1)]; 
        dt_phys_hist = [dt_phys_hist; zeros(grow, 1)]; 
        newton_it_hist = [newton_it_hist; zeros(grow, 1)]; %
        newton_res_hist = [newton_res_hist; zeros(grow, 1)]; %
    end

    x = x_next;
    Rk = x(1:d);
    Vk = x(d+1:2*d);
    tk = x(2*d+1);
    ptk = x(2*d+2);
    Hk = 0.5 * sum(Vk.^2) - mu / norm(Rk, 2);

    tau_hist(n) = tau_hist(n-1) + h_step;
    t_hist(n) = tk;
    R_hist(1, :, n) = Rk.';
    V_hist(1, :, n) = Vk.';
    pt_hist(n) = ptk;
    C_hist(n) = Hk + ptk;
    dt_phys_hist(n-1) = tk - t_hist(n-1);
    newton_it_hist(n-1) = nit;
    newton_res_hist(n-1) = fres;
end

tau_hist = tau_hist(1:n);
t_hist = t_hist(1:n);
R_hist = R_hist(:, :, 1:n);
V_hist = V_hist(:, :, 1:n);
pt_hist = pt_hist(1:n);
C_hist = C_hist(1:n);
if n > 1
    dt_phys_hist = dt_phys_hist(1:n-1);
    newton_it_hist = newton_it_hist(1:n-1);
    newton_res_hist = newton_res_hist(1:n-1);
else
    dt_phys_hist = zeros(0, 1);
    newton_it_hist = zeros(0, 1);
    newton_res_hist = zeros(0, 1);
end

stats = struct();
stats.accepted_steps = n - 1;
stats.h_tau = h_tau;
stats.rhs_evals = rhs_evals;
stats.pt_hist = pt_hist;
stats.constraint_hist = C_hist;
stats.dt_phys_hist = dt_phys_hist;
stats.newton_it_hist = newton_it_hist;
stats.newton_res_hist = newton_res_hist;
if ~isempty(dt_phys_hist)
    stats.dt_phys_min = min(dt_phys_hist);
    stats.dt_phys_max = max(dt_phys_hist);
else
    stats.dt_phys_min = NaN;
    stats.dt_phys_max = NaN;
end
if ~isempty(newton_it_hist)
    stats.newton_it_max = max(newton_it_hist);
    stats.newton_res_max = max(newton_res_hist);
else
    stats.newton_it_max = 0;
    stats.newton_res_max = 0.0;
end
end

function [x_next, nit, fres, ok, evals_step] = midpoint_step_newton( ...
    x, h, mu, g_fun, grad_g_fun, g_min, g_max, tol, maxit)
n = numel(x);
x_next = x;  % predictor
evals_step = 0;
ok = false;

for nit = 1:maxit
    [F, evals_f] = midpoint_residual(x_next, x, h, mu, g_fun, grad_g_fun, g_min, g_max);
    evals_step = evals_step + evals_f;
    fres = norm(F, inf);
    if fres < tol
        ok = true;
        return;
    end

    J = zeros(n, n);
    for j = 1:n
        dx = zeros(n, 1);
        eps_j = 1e-8 * max(1.0, abs(x_next(j)));
        dx(j) = eps_j;
        [Fp, evals_fp] = midpoint_residual(x_next + dx, x, h, mu, g_fun, grad_g_fun, g_min, g_max);
        evals_step = evals_step + evals_fp;
        J(:, j) = (Fp - F) / eps_j;
    end

    delta = -J \ F;
    x_next = x_next + delta;

    if norm(delta, inf) < tol * max(1.0, norm(x_next, inf))
        ok = true;
        [F2, evals_f2] = midpoint_residual(x_next, x, h, mu, g_fun, grad_g_fun, g_min, g_max);
        evals_step = evals_step + evals_f2;
        fres = norm(F2, inf);
        return;
    end
end

[F_end, evals_fe] = midpoint_residual(x_next, x, h, mu, g_fun, grad_g_fun, g_min, g_max);
evals_step = evals_step + evals_fe;
fres = norm(F_end, inf);
end

function [F, evals] = midpoint_residual(y, x, h, mu, g_fun, grad_g_fun, g_min, g_max)
xm = 0.5 * (x + y);
f = extended_rhs(xm, mu, g_fun, grad_g_fun, g_min, g_max);
F = y - x - h * f;
evals = 1;
end

function f = extended_rhs(x, mu, g_fun, grad_g_fun, g_min, g_max)
d = (numel(x) - 2) / 2;
R = x(1:d);
V = x(d+1:2*d);
pt = x(2*d+2);

r = norm(R, 2);
if r <= 0
    error('Encountered zero radius in canonical_sundman_midpoint.');
end

H = 0.5 * sum(V.^2) - mu / r;
g = clamp_scalar(g_fun(R), g_min, g_max);
gradg = grad_g_fun(R);
gradg = gradg(:);
if numel(gradg) ~= d
    error('grad_g_fun must return a vector of length %d.', d);
end

a_phys = -mu * R / r^3;
dR = g * V;
dV = g * a_phys - (H + pt) * gradg;
dt = g;
dpt = 0.0;

f = [dR; dV; dt; dpt];
end

function y = clamp_scalar(x, xmin, xmax)
y = min(xmax, max(xmin, x));
end

function v = get_opt(opts, name, default_value)
if isfield(opts, name)
    v = opts.(name);
else
    v = default_value;
end
end
