clear;
clc;

% Long-time comparison: velocity-Verlet (symplectic, order 2) vs BS-RK23.
% Metrics:
%   - orbit geometry over long integration
%   - energy drift DeltaE(t)
%   - computational cost (function evaluations)

params = struct('G', 1.0, 'm0', 1.0, 'masses', 1.0, 'dim', 2, 'eps_soft', 0.0);
params = nbody_prepare_params(params);

r0 = 1.0;
mu = params.G * params.m0;
v_circ = sqrt(mu / r0);
omega = sqrt(mu / r0^3);

% Recompute the part-2.1.1 reference step from the 1% circular-orbit target.
h_ref = find_reference_step(params, r0, omega, v_circ);
h_list = [h_ref, 1.5 * h_ref, 2 * h_ref];
n_periods = 100;

cases = struct([]);
cases(1).name = 'circular';
cases(1).v_factor = 1.0;
cases(1).has_exact = true;

cases(2).name = 'elliptical';
cases(2).v_factor = 0.9;
cases(2).has_exact = false;

results = struct([]);

fprintf('case\t\t h\t\t periods\t method\t\t max|DeltaE|\t\t cost_eval\t [r_min,r_max]\t max|r-1|\n');
for ic = 1:numel(cases)
    vy0 = cases(ic).v_factor * v_circ;
    T = orbital_period_from_initial(r0, vy0, mu);
    R0 = [r0, 0.0];
    V0 = [0.0, vy0];
    y0 = pack_state(R0, V0);
    tf_target = n_periods * T;

    for ih = 1:numel(h_list)
        h = h_list(ih);
        tf = tf_target;
        n_period_eff = n_periods;

        [tv, Rv, Vv, sv] = verlet_fixed(@nbody_acceleration, [0, tf], R0, V0, h, params);
        Ev = energy_series_from_rv(Rv, Vv, params);
        dEv = Ev - Ev(1);
        max_dEv = max(abs(dEv));

        [tb, Yb, sb] = bs23_fixed(@nbody_rhs, [0, tf], y0, h, params);
        [Rb, Vb] = rv_from_state_history(Yb, params);
        Eb = energy_series_from_rv(Rb, Vb, params);
        dEb = Eb - Eb(1);
        max_dEb = max(abs(dEb));

        rv = squeeze(Rv(1, :, :));
        rb = squeeze(Rb(1, :, :));
        rnorm_v = sqrt(sum(rv.^2, 1));
        rnorm_b = sqrt(sum(rb.^2, 1));
        rdev_v = max(abs(rnorm_v - r0));
        rdev_b = max(abs(rnorm_b - r0));

        fprintf('%-9s\t %.5f\t %.1f\t Verlet\t\t %.3e\t %.0f\t\t [%.3f, %.3f]\t %.3e\n', ...
            cases(ic).name, h, n_period_eff, max_dEv, sv.accel_evals, min(rnorm_v), max(rnorm_v), rdev_v);
        fprintf('%-9s\t %.5f\t %.1f\t BS23\t\t %.3e\t %.0f\t\t [%.3f, %.3f]\t %.3e\n', ...
            cases(ic).name, h, n_period_eff, max_dEb, sb.nfev, min(rnorm_b), max(rnorm_b), rdev_b);

        results(ic, ih).case_name = cases(ic).name; %
        results(ic, ih).h = h;
        results(ic, ih).T = T;
        results(ic, ih).n_period_eff = n_period_eff;
        results(ic, ih).tv = tv;
        results(ic, ih).tb = tb;
        results(ic, ih).Rv = Rv;
        results(ic, ih).Rb = Rb;
        results(ic, ih).dEv = dEv;
        results(ic, ih).dEb = dEb;
        results(ic, ih).max_dEv = max_dEv;
        results(ic, ih).max_dEb = max_dEb;
        results(ic, ih).cost_verlet = sv.accel_evals;
        results(ic, ih).cost_bs23 = sb.nfev;
        results(ic, ih).rdev_v = rdev_v;
        results(ic, ih).rdev_b = rdev_b;
    end
end

for ic = 1:numel(cases)
    % Plot orbit geometry for reference step size.
    ref = results(ic, 1);
    figure('Name', sprintf('%s orbit long-time (h=%.4f)', cases(ic).name, ref.h));
    hold on;
    traj_v = squeeze(ref.Rv(1, :, :));
    traj_b = squeeze(ref.Rb(1, :, :));
    plot(traj_v(1, :), traj_v(2, :), 'LineWidth', 1.2, 'DisplayName', 'Verlet');
    plot(traj_b(1, :), traj_b(2, :), 'LineWidth', 1.2, 'DisplayName', 'BS-RK23');
    if cases(ic).has_exact
        t_ex = linspace(0, ref.tv(end), 2000);
        plot(r0 * cos(omega * t_ex), r0 * sin(omega * t_ex), 'k--', ...
            'LineWidth', 1.0, 'DisplayName', 'exact');
    end
    plot(0, 0, 'kp', 'MarkerFaceColor', 'y', 'DisplayName', 'star');
    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');
    title(sprintf('%s orbit over %.1f periods (h=%.4f)', cases(ic).name, ref.n_period_eff, ref.h));
    legend('Location', 'best');

    % Plot long-time energy drift for reference step size.
    figure('Name', sprintf('%s energy drift (h=%.4f)', cases(ic).name, ref.h));
    plot(ref.tv, ref.dEv, 'LineWidth', 1.2, 'DisplayName', 'Verlet'); hold on;
    plot(ref.tb, ref.dEb, 'LineWidth', 1.2, 'DisplayName', 'BS-RK23');
    grid on;
    xlabel('t');
    ylabel('\Delta E(t)');
    title(sprintf('%s: energy drift over %.1f periods', cases(ic).name, ref.n_period_eff));
    legend('Location', 'best');

    % Robustness with larger h: max |DeltaE| versus h.
    max_de_v = zeros(1, numel(h_list));
    max_de_b = zeros(1, numel(h_list));
    for ih = 1:numel(h_list)
        max_de_v(ih) = results(ic, ih).max_dEv;
        max_de_b(ih) = results(ic, ih).max_dEb;
    end
    figure('Name', sprintf('%s robustness: max |DeltaE| vs h', cases(ic).name));
    loglog(h_list, max_de_v, 'o-', 'LineWidth', 1.4, 'DisplayName', 'Verlet'); hold on;
    loglog(h_list, max_de_b, 's-', 'LineWidth', 1.4, 'DisplayName', 'BS-RK23');
    grid on;
    xlabel('h');
    ylabel('max |\Delta E(t)|');
    title(sprintf('%s robustness over %.1f periods', cases(ic).name, n_periods));
    legend('Location', 'best');
end

function T = orbital_period_from_initial(r0, vy0, mu)
E = 0.5 * vy0^2 - mu / r0;
if E >= 0
    error('Initial condition is not a bound orbit (E >= 0).');
end
a = -mu / (2 * E);
T = 2 * pi * sqrt(a^3 / mu);
end

function h_ref = find_reference_step(params, r0, omega, v_circ)
T = 2 * pi / omega;
y0 = pack_state([r0, 0.0], [0.0, v_circ]);
target_err = 1e-2 * r0;
M_lo = 4;
M_hi = 16;

while true
    h = T / M_hi;
    [t, Y, ~] = bs23_fixed(@nbody_rhs, [0, T], y0, h, params);
    e = max_position_error_vs_exact(t, Y, r0, omega, params);
    if e <= target_err
        break;
    end
    M_hi = 2 * M_hi;
end

while M_lo + 1 < M_hi
    M_mid = floor((M_lo + M_hi) / 2);
    h = T / M_mid;
    [t, Y, ~] = bs23_fixed(@nbody_rhs, [0, T], y0, h, params);
    e = max_position_error_vs_exact(t, Y, r0, omega, params);
    if e <= target_err
        M_hi = M_mid;
    else
        M_lo = M_mid;
    end
end

h_ref = T / M_hi;
end

function err_max = max_position_error_vs_exact(t, Y, r0, omega, params)
R = reshape(Y(1:params.N*params.dim, :), params.N, params.dim, []);
R = squeeze(R(1, :, :));
R_ex = [r0 * cos(omega * t).'; r0 * sin(omega * t).'];
err_max = max(vecnorm(R - R_ex, 2, 1));
end

function E = energy_series_from_rv(R_hist, V_hist, params)
nt = size(R_hist, 3);
E = zeros(nt, 1);
for k = 1:nt
    E(k) = nbody_energy(pack_state(R_hist(:, :, k), V_hist(:, :, k)), params);
end
end

function [R_hist, V_hist] = rv_from_state_history(Y, params)
nt = size(Y, 2);
R_hist = zeros(params.N, params.dim, nt);
V_hist = zeros(params.N, params.dim, nt);
for k = 1:nt
    [Rk, Vk] = unpack_state(Y(:, k), params.N, params.dim);
    R_hist(:, :, k) = Rk;
    V_hist(:, :, k) = Vk;
end
end
