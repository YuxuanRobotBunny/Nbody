clear;
clc;

% Project A:
% Long-time energy behavior comparison:
%   1) fixed-step Verlet (symplectic)
%   2) fixed-step BS-RK23 (non-symplectic)
%   3) adaptive BS-RK23 (variable-step, non-symplectic)
%   4) adaptive naive Verlet (variable-step, generally non-symplectic)

params = struct('G', 1.0, 'm0', 1.0, 'masses', 1.0, 'dim', 2, 'eps_soft', 0.0);
params = nbody_prepare_params(params);

r0 = 1.0;
mu = params.G * params.m0;
v_circ = sqrt(mu / r0);
v_factor = 0.5;
vy0 = v_factor * v_circ;

T = orbital_period_from_initial(r0, vy0, mu);
n_periods = 200;
t_end = n_periods * T;

R0 = [r0, 0.0];
V0 = [0.0, vy0];
y0 = pack_state(R0, V0);

h = T / 300;
h_cost = 4 * h;   % approximately matches cost with Verlet at h

fprintf('Setup: v(0)=%.3f v_c, T=%.6f, t_end=%.6f (%.1f periods)\n', ...
    v_factor, T, t_end, n_periods);

% Fixed-step Verlet (symplectic).
[tv, Rv, Vv, sv] = verlet_fixed(@nbody_acceleration, [0, t_end], R0, V0, h, params);
[dEv, rV] = diagnostics_from_rv(tv, Rv, Vv, params);

% Fixed-step BS23, same h.
[tb_same, Yb_same, sb_same] = bs23_fixed(@nbody_rhs, [0, t_end], y0, h, params);
[Rb_same, Vb_same] = rv_from_state_history(Yb_same, params);
[dEb_same, rB_same] = diagnostics_from_rv(tb_same, Rb_same, Vb_same, params);

% Fixed-step BS23, approximately cost-matched (h multiplied by ~4).
[tb_cost, Yb_cost, sb_cost] = bs23_fixed(@nbody_rhs, [0, t_end], y0, h_cost, params);
[Rb_cost, Vb_cost] = rv_from_state_history(Yb_cost, params);
[dEb_cost, rB_cost] = diagnostics_from_rv(tb_cost, Rb_cost, Vb_cost, params);

% Adaptive BS23 (variable-step, improved controller).
tol = 1e-4;
opts_bs = struct('err_mode', 'full', 'rtol', tol, 'atol', 1e-10, ...
    'controller', 'pi', 'fac_max', 2.0, 'fac_max_reject', 0.8, ...
    'h_min', 1e-8, 'h_max', 0.2);
[tb_ad, Yb_ad, sb_ad] = bs23_adaptive(@nbody_rhs, [0, t_end], y0, 0.02, tol, params, opts_bs);
[Rb_ad, Vb_ad] = rv_from_state_history(Yb_ad, params);
[dEb_ad, rB_ad] = diagnostics_from_rv(tb_ad, Rb_ad, Vb_ad, params);

% Adaptive naive Verlet (variable-step, intentionally non-symplectic).
opts_vad = struct('rtol', tol, 'atol', 1e-10, 'fac_max', 2.0, 'fac_max_reject', 0.8, ...
    'h_min', 1e-8, 'h_max', 0.2);
[tv_ad, Rv_ad, Vv_ad, sv_ad] = verlet_adaptive_naive(@nbody_acceleration, ...
    [0, t_end], R0, V0, 0.02, tol, params, opts_vad);
[dEv_ad, rV_ad] = diagnostics_from_rv(tv_ad, Rv_ad, Vv_ad, params);

fprintf('\nMethod summary:\n');
print_summary('Verlet fixed (h)', max(abs(dEv)), min(rV), max(rV), sv.accel_evals);
print_summary('BS23 fixed (h)', max(abs(dEb_same)), min(rB_same), max(rB_same), sb_same.nfev);
print_summary('BS23 fixed (4h, cost-match approx)', max(abs(dEb_cost)), min(rB_cost), max(rB_cost), sb_cost.nfev);
print_summary('BS23 adaptive', max(abs(dEb_ad)), min(rB_ad), max(rB_ad), sb_ad.nfev);
print_summary('Verlet adaptive naive', max(abs(dEv_ad)), min(rV_ad), max(rV_ad), sv_ad.accel_evals);

fprintf('\nCost ratio (relative to Verlet fixed):\n');
fprintf('BS23 fixed (h): %.2fx\n', sb_same.nfev / sv.accel_evals);
fprintf('BS23 fixed (4h): %.2fx\n', sb_cost.nfev / sv.accel_evals);
fprintf('BS23 adaptive: %.2fx\n', sb_ad.nfev / sv.accel_evals);
fprintf('Verlet adaptive naive: %.2fx\n', sv_ad.accel_evals / sv.accel_evals);

% -------------------------------------------------------------------------
% Plot 1: DeltaE(t)
% -------------------------------------------------------------------------
figure('Name', 'Project A: DeltaE(t) comparison');
plot(tv, dEv, 'LineWidth', 1.5, 'DisplayName', 'Verlet fixed (h)'); hold on;
plot(tb_same, dEb_same, 'LineWidth', 1.2, 'DisplayName', 'BS23 fixed (h)');
plot(tb_cost, dEb_cost, 'LineWidth', 1.2, 'DisplayName', 'BS23 fixed (4h, cost-match approx)');
plot(tb_ad, dEb_ad, 'LineWidth', 1.2, 'DisplayName', 'BS23 adaptive');
plot(tv_ad, dEv_ad, 'LineWidth', 1.2, 'DisplayName', 'Verlet adaptive naive');
grid on;
xlabel('t');
ylabel('\Delta E(t) = E(t)-E(0)');
title(sprintf('Energy error over %.1f periods (v(0)=%.3f v_c)', n_periods, v_factor));
legend('Location', 'best');

% -------------------------------------------------------------------------
% Plot 2: Radius history r(t)
% -------------------------------------------------------------------------
figure('Name', 'Project A: radius drift comparison');
plot(tv, rV, 'LineWidth', 1.5, 'DisplayName', 'Verlet fixed (h)'); hold on;
plot(tb_same, rB_same, 'LineWidth', 1.2, 'DisplayName', 'BS23 fixed (h)');
plot(tb_cost, rB_cost, 'LineWidth', 1.2, 'DisplayName', 'BS23 fixed (4h, cost-match approx)');
plot(tb_ad, rB_ad, 'LineWidth', 1.2, 'DisplayName', 'BS23 adaptive');
plot(tv_ad, rV_ad, 'LineWidth', 1.2, 'DisplayName', 'Verlet adaptive naive');
grid on;
xlabel('t');
ylabel('r(t)=||R(t)||_2');
title('Orbital radius drift comparison');
legend('Location', 'best');

% -------------------------------------------------------------------------
% Plot 3: Accepted step-size history for adaptive methods
% -------------------------------------------------------------------------
figure('Name', 'Project A: adaptive step-size history');
acc_b = sb_ad.accepted_attempt;
acc_v = sv_ad.accepted_attempt;
plot(sb_ad.t_attempt(acc_b), sb_ad.h_attempt(acc_b), '.-', 'DisplayName', 'BS23 adaptive accepted'); hold on;
plot(sv_ad.t_attempt(acc_v), sv_ad.h_attempt(acc_v), '.-', 'DisplayName', 'Verlet adaptive naive accepted');
grid on;
xlabel('t');
ylabel('accepted h');
title('Variable step-size history');
legend('Location', 'best');

function T = orbital_period_from_initial(r0, vy0, mu)
E = 0.5 * vy0^2 - mu / r0;
if E >= 0
    error('Initial condition is not a bound orbit (E >= 0).');
end
a = -mu / (2 * E);
T = 2 * pi * sqrt(a^3 / mu);
end

function [dE, rnorm] = diagnostics_from_rv(t, R_hist, V_hist, params)
nt = numel(t);
E = zeros(nt, 1);
rnorm = zeros(nt, 1);
for k = 1:nt
    Rk = R_hist(:, :, k);
    Vk = V_hist(:, :, k);
    E(k) = nbody_energy(pack_state(Rk, Vk), params);
    rnorm(k) = norm(Rk(1, :), 2);
end
dE = E - E(1);
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

function print_summary(name, max_dE, rmin, rmax, n_eval)
fprintf('%-34s  max|DeltaE|=%-11.3e  r in [%.6f, %.6f]  eval=%d\n', ...
    name, max_dE, rmin, rmax, n_eval);
end
