clear;
clc;

% This script is a teaching-oriented view of adaptive BS-RK23 behavior.
% It explains:
%   (1) what x(t) curves mean,
%   (2) how adaptive step size reacts to orbital difficulty.

params = struct('G', 1.0, 'm0', 1.0, 'masses', 1.0, 'dim', 2, 'eps_soft', 0.0);
params = nbody_prepare_params(params);

r0 = 1.0;
mu = params.G * params.m0;
v_circ = sqrt(mu / r0);

velocity_factors = [1.0, 0.5, 0.25, 0.125];
tol = 1e-6;
h0 = 0.02;
opts = struct( ...
    'err_mode', 'full', ...
    'rtol', tol, ...
    'atol', 1e-10, ...
    'controller', 'pi', ...
    'fac_max', 2.0, ...
    'fac_max_reject', 0.8);

results = struct([]);

fprintf('factor\tT\t\tmin(r)\t\tmax(r)\t\tmax(|x|)\tacc\trej\trej_rate\tmin(h_acc)\tmax(h_acc)\n');
for i = 1:numel(velocity_factors)
    f = velocity_factors(i);
    vy0 = f * v_circ;
    T = orbital_period_from_initial(r0, vy0, mu);
    y0 = pack_state([r0, 0.0], [0.0, vy0]);

    [t, Y, stats] = bs23_adaptive(@nbody_rhs, [0, T], y0, h0, tol, params, opts);

    R = reshape(Y(1:params.N*params.dim, :), params.N, params.dim, []);
    x = squeeze(R(1, 1, :));
    y = squeeze(R(1, 2, :));
    r = sqrt(x.^2 + y.^2);

    h_acc = stats.h_attempt(stats.accepted_attempt);
    rej_rate = stats.rejected_steps / max(1, (stats.accepted_steps + stats.rejected_steps));

    fprintf('%.3f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\t%.3f\t\t%.3e\t%.3e\n', ...
        f, T, min(r), max(r), max(abs(x)), stats.accepted_steps, stats.rejected_steps, ...
        rej_rate, min(h_acc), max(h_acc));

    results(i).f = f; %
    results(i).T = T;
    results(i).t = t;
    results(i).x = x;
    results(i).y = y;
    results(i).r = r;
    results(i).stats = stats;
end

% -------------------------------------------------------------------------
% Figure 1: Circular case only (numerical vs exact x(t))
% -------------------------------------------------------------------------
res_circ = results(1);  % f=1.0
t_circ = res_circ.t;
x_num = res_circ.x;
x_exact = r0 * cos(t_circ);

figure('Name', 'Adaptive BS-RK23 Explained: Circular Orbit');
plot(t_circ, x_num, 'b-', 'LineWidth', 1.5, 'DisplayName', 'numerical x(t)'); hold on;
plot(t_circ, x_exact, 'k--', 'LineWidth', 1.3, 'DisplayName', 'exact x(t)=cos(t)');
grid on;
xlabel('t');
ylabel('x(t)');
title(sprintf('Circular Orbit, rtol=%.0e, atol=%.0e: numerical vs exact', opts.rtol, opts.atol));
legend('Location', 'best');

% -------------------------------------------------------------------------
% Figure 2: Compare all velocity factors on normalized time tau=t/T
% -------------------------------------------------------------------------
figure('Name', 'Adaptive BS-RK23 Explained: x(t) by Eccentricity');
hold on;
for i = 1:numel(results)
    tau = results(i).t / results(i).T;
    plot(tau, results(i).x, 'LineWidth', 1.3, ...
        'DisplayName', sprintf('v(0)=%.3f v_c', results(i).f));
end
grid on;
xlabel('\tau = t/T (one orbit normalized to [0,1])');
ylabel('x(t)');
title(sprintf('x(t) Comparison (rtol=%.0e, atol=%.0e)', opts.rtol, opts.atol));
legend('Location', 'best');

% -------------------------------------------------------------------------
% Figure 3: Link orbit difficulty to adaptive step size for hardest case
% -------------------------------------------------------------------------
res_hard = results(end);  % f=0.125
[~, idx_peri] = min(res_hard.r);
t_peri = res_hard.t(idx_peri);
tau_orbit = res_hard.t / res_hard.T;
tau_attempt = res_hard.stats.t_attempt / res_hard.T;

acc = res_hard.stats.accepted_attempt;
rej = ~acc;

figure('Name', 'Adaptive BS-RK23 Explained: r(t) and h(t)');
subplot(2, 1, 1);
plot(tau_orbit, res_hard.r, 'LineWidth', 1.4);
grid on;
xlabel('\tau = t/T');
ylabel('r(t)');
title('Hard case (v(0)=0.125 v_c): radius over one orbit');
xline(t_peri / res_hard.T, 'r--', 'periapsis', 'LineWidth', 1.2, 'LabelVerticalAlignment', 'middle');

subplot(2, 1, 2);
plot(tau_attempt(acc), res_hard.stats.h_attempt(acc), 'b.-', 'DisplayName', 'accepted'); hold on;
plot(tau_attempt(rej), res_hard.stats.h_attempt(rej), 'rx', 'DisplayName', 'rejected');
grid on;
xlabel('\tau = t/T');
ylabel('attempted h');
title('Adaptive step-size history (small h near periapsis is expected)');
legend('Location', 'best');
xline(t_peri / res_hard.T, 'r--', 'periapsis', 'LineWidth', 1.2, 'LabelVerticalAlignment', 'middle');

% -------------------------------------------------------------------------
% Quick interpretation text in command window
% -------------------------------------------------------------------------
fprintf('\nInterpretation guide:\n');
fprintf('1) Circular case should track cos(t) closely.\n');
fprintf('2) Smaller v(0) => more eccentric orbit => stronger acceleration near periapsis.\n');
fprintf('3) Near periapsis, adaptive solver should reduce h and may reject steps.\n');
fprintf('4) If x(t) amplitude becomes physically implausible, tolerance is too loose for that case.\n');

function T = orbital_period_from_initial(r0, vy0, mu)
E = 0.5 * vy0^2 - mu / r0;
if E >= 0
    error('Initial condition is not a bound orbit (E >= 0).');
end
a = -mu / (2 * E);
T = 2 * pi * sqrt(a^3 / mu);
end
