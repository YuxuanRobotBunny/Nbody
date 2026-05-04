function [t_hist, uhat_hist, stats] = etdrk2_kdv(grid, tspan, u0_hat, dt)
%ETDRK2_KDV Second-order exponential RK integrator for KdV.

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
while t_hist(n) < tf - t_eps
    dt_step = min(dt, tf - t_hist(n));
    [E, phi1, phi2] = etdrk2_coeffs(grid.lin, dt_step);

    uhat_n = uhat_hist(:, n);
    Bn = kdv_nonlinear_hat(uhat_n, grid);
    uhat_star = E .* uhat_n + dt_step * phi1 .* Bn;
    Bstar = kdv_nonlinear_hat(uhat_star, grid);
    uhat_next = uhat_star + dt_step * phi2 .* (Bstar - Bn);
    nonlinear_evals = nonlinear_evals + 2;

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
    'method', 'ETDRK2');
end

function [E, phi1, phi2] = etdrk2_coeffs(L, dt)
z = dt * L;
E = exp(z);
phi1 = phi1_stable(z);
phi2 = phi2_stable(z);
end

function y = phi1_stable(z)
y = zeros(size(z));
small = abs(z) < 1e-7;
zs = z(small);
y(small) = 1 + zs / 2 + zs.^2 / 6 + zs.^3 / 24;
y(~small) = (exp(z(~small)) - 1) ./ z(~small);
end

function y = phi2_stable(z)
y = zeros(size(z));
small = abs(z) < 1e-5;
zs = z(small);
y(small) = 0.5 + zs / 6 + zs.^2 / 24 + zs.^3 / 120;
y(~small) = (exp(z(~small)) - 1 - z(~small)) ./ (z(~small).^2);
end
