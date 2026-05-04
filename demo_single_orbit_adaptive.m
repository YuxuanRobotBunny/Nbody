clear;
clc;

% Single-planet elliptical orbits: adaptive BS-RK23 validation.
params = struct('G', 1.0, 'm0', 1.0, 'masses', 1.0, 'dim', 2, 'eps_soft', 0.0);
params = nbody_prepare_params(params);

r0 = 1.0;
mu = params.G * params.m0;
v_circ = sqrt(mu / r0);
target_abs_err = 1e-2 * r0;

velocity_factors = [1.0, 0.5, 0.25, 0.125];
h0 = 0.02;

fprintf('factor\t atol\t\t final_pos_err_inf\t target_met\t accepted\t rejected\t nfev\n');
for i = 1:numel(velocity_factors)
    f = velocity_factors(i);
    vy0 = f * v_circ;
    T = orbital_period_from_initial(r0, vy0, mu);
    y0 = pack_state([r0, 0.0], [0.0, vy0]);

    [atol_star, final_err, stats] = find_position_atol_for_target(y0, T, h0, target_abs_err, params);
    target_met = final_err <= target_abs_err;
    fprintf('%.3f\t %.1e\t %.3e\t\t %d\t\t %d\t\t %d\t\t %d\n', ...
        f, atol_star, final_err, target_met, stats.accepted_steps, stats.rejected_steps, stats.nfev);
end

% Plot x(t) for different eccentricities using the calibrated absolute tolerance.
figure;
hold on;
for i = 1:numel(velocity_factors)
    f = velocity_factors(i);
    vy0 = f * v_circ;
    T = orbital_period_from_initial(r0, vy0, mu);
    y0 = pack_state([r0, 0.0], [0.0, vy0]);
    [atol_star, ~, ~] = find_position_atol_for_target(y0, T, h0, target_abs_err, params);
    opts = adaptive_position_opts(atol_star);
    [t, Y, ~] = bs23_adaptive(@nbody_rhs, [0, T], y0, h0, atol_star, params, opts);

    R = reshape(Y(1:params.N*params.dim, :), params.N, params.dim, []);
    x = squeeze(R(1, 1, :));
    plot(t, x, 'LineWidth', 1.2, 'DisplayName', sprintf('v(0)=%.3f v_c, atol=%.1e', f, atol_star));
end
grid on;
xlabel('t');
ylabel('x(t)');
title(sprintf('Adaptive BS-RK23: x(t), target final |error| <= %.0e', target_abs_err));
legend('Location', 'best');

% Plot attempted step sizes including rejected steps for most eccentric orbit.
f = velocity_factors(end);
vy0 = f * v_circ;
T = orbital_period_from_initial(r0, vy0, mu);
y0 = pack_state([r0, 0.0], [0.0, vy0]);
[atol_star, ~, ~] = find_position_atol_for_target(y0, T, h0, target_abs_err, params);
opts = adaptive_position_opts(atol_star);
[~, ~, stats] = bs23_adaptive(@nbody_rhs, [0, T], y0, h0, atol_star, params, opts);

acc = stats.accepted_attempt;
rej = ~acc;

figure;
plot(stats.t_attempt(acc), stats.h_attempt(acc), 'b.-', 'DisplayName', 'accepted'); hold on;
plot(stats.t_attempt(rej), stats.h_attempt(rej), 'rx', 'DisplayName', 'rejected');
grid on;
xlabel('t');
ylabel('attempted h');
title(sprintf('Adaptive Step Sizes (v(0)=%.3f v_c, atol=%.1e)', f, atol_star));
legend('Location', 'best');

function T = orbital_period_from_initial(r0, vy0, mu)
E = 0.5 * vy0^2 - mu / r0;
if E >= 0
    error('Initial condition is not a bound orbit (E >= 0).');
end
a = -mu / (2 * E);
T = 2 * pi * sqrt(a^3 / mu);
end

function opts = adaptive_position_opts(atol_pos)
opts = struct('err_mode', 'position', 'rtol', 0.0, 'atol', atol_pos, ...
    'controller', 'pi', 'fac_max', 2.0, 'fac_max_reject', 0.8);
end

function [atol_star, final_err, stats_star] = find_position_atol_for_target(y0, T, h0, target_abs_err, params)
atol_hi = target_abs_err;
[final_err_hi, stats_hi] = run_with_atol(y0, T, h0, atol_hi, params);
if final_err_hi <= target_abs_err
    atol_star = atol_hi;
    final_err = final_err_hi;
    stats_star = stats_hi;
    return;
end

atol_lo = atol_hi;
stats_lo = stats_hi;
final_err_lo = final_err_hi;
while final_err_lo > target_abs_err
    atol_lo = 0.5 * atol_lo;
    [final_err_lo, stats_lo] = run_with_atol(y0, T, h0, atol_lo, params);
    if atol_lo < 1e-12
        error('Failed to meet final absolute target %.3e.', target_abs_err);
    end
end

for it = 1:18
    atol_mid = sqrt(atol_lo * atol_hi);
    [final_err_mid, stats_mid] = run_with_atol(y0, T, h0, atol_mid, params);
    if final_err_mid <= target_abs_err
        atol_lo = atol_mid;
        final_err_lo = final_err_mid;
        stats_lo = stats_mid;
    else
        atol_hi = atol_mid;
    end
end

atol_star = atol_lo;
final_err = final_err_lo;
stats_star = stats_lo;
end

function [final_err, stats] = run_with_atol(y0, T, h0, atol_pos, params)
opts = adaptive_position_opts(atol_pos);
[~, Y, stats] = bs23_adaptive(@nbody_rhs, [0, T], y0, h0, atol_pos, params, opts);
[R_end, ~] = unpack_state(Y(:, end), params.N, params.dim);
final_err = norm(R_end(1, :) - [1.0, 0.0], inf);
end
