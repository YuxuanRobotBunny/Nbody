clear;
clc;

% Project C:
% Phase-space area test for harmonic oscillator q' = p, p' = -omega^2 q.
%
% We compare:
%   1) fixed-step Verlet (symplectic)
%   2) fixed-step RK23 (non-symplectic)
%   3) state-dependent-step RK23 (non-symplectic)
%   4) state-dependent-step Verlet (generally non-symplectic)
%
% The state-dependent step size is chosen as:
%   h_i = h0 / (1 + gamma*(q_i^2 + p_i^2)),
% with clipping. Because h depends on state, this map is not a standard
% symplectic one-step map.

omega = 1.0;
h0 = 0.08;
h_min = 0.01;
h_max = 0.12;
gamma = 10.0;
n_steps = 5000;
record_stride = 10;

% Cloud in (q,p): choose a nontrivial patch so state-dependent h varies visibly.
q_c = 1.0;
p_c = 0.0;
dq = 0.18;
dp = 0.18;
n_cloud = 260;
theta = linspace(0, 2*pi, n_cloud + 1);
theta(end) = [];

Q0 = q_c + dq * cos(theta);
P0 = p_c + dp * sin(theta);
A0 = cloud_area(Q0, P0);

[area_v, t_v, Qv, Pv] = track_area_verlet_fixed(Q0, P0, omega, h0, n_steps, record_stride);
[area_rk, t_rk, Qrk, Prk] = track_area_rk23_fixed(Q0, P0, omega, h0, n_steps, record_stride);
[area_rk_sd, t_rk_sd, Qrk_sd, Prk_sd, hmean_rk_sd] = ...
    track_area_rk23_state_dep(Q0, P0, omega, h0, n_steps, record_stride, gamma, h_min, h_max);
[area_v_sd, t_v_sd, Qv_sd, Pv_sd, hmean_v_sd] = ...
    track_area_verlet_state_dep(Q0, P0, omega, h0, n_steps, record_stride, gamma, h_min, h_max);

fprintf('Initial cloud area A0 = %.10e\n', A0);
fprintf('Final area ratios A/A0:\n');
fprintf('  Verlet fixed                : %.8f\n', area_v(end) / A0);
fprintf('  RK23 fixed                  : %.8f\n', area_rk(end) / A0);
fprintf('  RK23 state-dependent h      : %.8f\n', area_rk_sd(end) / A0);
fprintf('  Verlet state-dependent h    : %.8f\n', area_v_sd(end) / A0);

% -------------------------------------------------------------------------
% Plot 1: area ratio
% -------------------------------------------------------------------------
figure('Name', 'Project C: phase-space area ratio');
plot(t_v, area_v / A0, 'LineWidth', 1.5, 'DisplayName', 'Verlet fixed'); hold on;
plot(t_rk, area_rk / A0, 'LineWidth', 1.2, 'DisplayName', 'RK23 fixed');
plot(t_rk_sd, area_rk_sd / A0, 'LineWidth', 1.2, 'DisplayName', 'RK23 state-dependent h');
plot(t_v_sd, area_v_sd / A0, 'LineWidth', 1.2, 'DisplayName', 'Verlet state-dependent h');
yline(1.0, 'k--', 'A/A_0 = 1');
grid on;
xlabel('pseudo-time t = n h_0');
ylabel('A(t)/A(0)');
title('Phase-space area preservation test');
legend('Location', 'best');

% -------------------------------------------------------------------------
% Plot 2: mean step size in state-dependent methods
% -------------------------------------------------------------------------
figure('Name', 'Project C: mean adaptive step');
plot(t_rk_sd, hmean_rk_sd, 'LineWidth', 1.2, 'DisplayName', 'RK23 mean h'); hold on;
plot(t_v_sd, hmean_v_sd, 'LineWidth', 1.2, 'DisplayName', 'Verlet mean h');
grid on;
xlabel('pseudo-time t = n h_0');
ylabel('mean h over cloud');
title('State-dependent step-size evolution');
legend('Location', 'best');

% -------------------------------------------------------------------------
% Plot 3: cloud snapshots
% -------------------------------------------------------------------------
figure('Name', 'Project C: cloud snapshots');
idx_v = [1, floor(numel(t_v)/2), numel(t_v)];
idx_vsd = [1, floor(numel(t_v_sd)/2), numel(t_v_sd)];

subplot(1, 2, 1);
hold on;
for j = 1:numel(idx_v)
    q = Qv(:, idx_v(j)); p = Pv(:, idx_v(j));
    k = convhull(q, p);
    plot(q(k), p(k), 'LineWidth', 1.2, 'DisplayName', sprintf('t=%.1f', t_v(idx_v(j))));
end
grid on;
xlabel('q'); ylabel('p');
title('Verlet fixed');
legend('Location', 'best');
axis equal;

subplot(1, 2, 2);
hold on;
for j = 1:numel(idx_vsd)
    q = Qv_sd(:, idx_vsd(j)); p = Pv_sd(:, idx_vsd(j));
    k = convhull(q, p);
    plot(q(k), p(k), 'LineWidth', 1.2, 'DisplayName', sprintf('t=%.1f', t_v_sd(idx_vsd(j))));
end
grid on;
xlabel('q'); ylabel('p');
title('Verlet state-dependent h');
legend('Location', 'best');
axis equal;

% ------------------------------ helpers ----------------------------------
function [A_hist, t_hist, Q_hist, P_hist] = track_area_verlet_fixed(Q0, P0, omega, h, n_steps, stride)
Q = Q0(:).';
P = P0(:).';
[A_hist, t_hist, Q_hist, P_hist] = init_record(Q, P, n_steps, stride);
rec = 1;
t = 0;
for n = 1:n_steps
    [Q, P] = verlet_step(Q, P, h, omega);
    t = t + h;
    if mod(n, stride) == 0
        rec = rec + 1;
        [A_hist, t_hist, Q_hist, P_hist] = save_record(A_hist, t_hist, Q_hist, P_hist, rec, t, Q, P);
    end
end
end

function [A_hist, t_hist, Q_hist, P_hist] = track_area_rk23_fixed(Q0, P0, omega, h, n_steps, stride)
Q = Q0(:).';
P = P0(:).';
[A_hist, t_hist, Q_hist, P_hist] = init_record(Q, P, n_steps, stride);
rec = 1;
t = 0;
for n = 1:n_steps
    [Q, P] = rk23_step(Q, P, h, omega);
    t = t + h;
    if mod(n, stride) == 0
        rec = rec + 1;
        [A_hist, t_hist, Q_hist, P_hist] = save_record(A_hist, t_hist, Q_hist, P_hist, rec, t, Q, P);
    end
end
end

function [A_hist, t_hist, Q_hist, P_hist, hmean_hist] = ...
    track_area_rk23_state_dep(Q0, P0, omega, h0, n_steps, stride, gamma, h_min, h_max)
Q = Q0(:).';
P = P0(:).';
[A_hist, t_hist, Q_hist, P_hist] = init_record(Q, P, n_steps, stride);
hmean_hist = zeros(size(A_hist));
hmean_hist(1) = h0;
rec = 1;
t = 0;
for n = 1:n_steps
    H = state_dependent_h(Q, P, h0, gamma, h_min, h_max);
    [Q, P] = rk23_step(Q, P, H, omega);
    t = t + h0;
    if mod(n, stride) == 0
        rec = rec + 1;
        [A_hist, t_hist, Q_hist, P_hist] = save_record(A_hist, t_hist, Q_hist, P_hist, rec, t, Q, P);
        hmean_hist(rec) = mean(H);
    end
end
hmean_hist = hmean_hist(1:rec);
end

function [A_hist, t_hist, Q_hist, P_hist, hmean_hist] = ...
    track_area_verlet_state_dep(Q0, P0, omega, h0, n_steps, stride, gamma, h_min, h_max)
Q = Q0(:).';
P = P0(:).';
[A_hist, t_hist, Q_hist, P_hist] = init_record(Q, P, n_steps, stride);
hmean_hist = zeros(size(A_hist));
hmean_hist(1) = h0;
rec = 1;
t = 0;
for n = 1:n_steps
    H = state_dependent_h(Q, P, h0, gamma, h_min, h_max);
    [Q, P] = verlet_step(Q, P, H, omega);
    t = t + h0;
    if mod(n, stride) == 0
        rec = rec + 1;
        [A_hist, t_hist, Q_hist, P_hist] = save_record(A_hist, t_hist, Q_hist, P_hist, rec, t, Q, P);
        hmean_hist(rec) = mean(H);
    end
end
hmean_hist = hmean_hist(1:rec);
end

function H = state_dependent_h(Q, P, h0, gamma, h_min, h_max)
H = h0 ./ (1 + gamma * (Q.^2 + P.^2));
H = min(h_max, max(h_min, H));
end

function [Qn, Pn] = verlet_step(Q, P, h, omega)
P_half = P - 0.5 .* h .* omega^2 .* Q;
Qn = Q + h .* P_half;
Pn = P_half - 0.5 .* h .* omega^2 .* Qn;
end

function [Q3, P3] = rk23_step(Q, P, h, omega)
k1q = P;
k1p = -omega^2 .* Q;

q2 = Q + h .* (0.5 .* k1q);
p2 = P + h .* (0.5 .* k1p);
k2q = p2;
k2p = -omega^2 .* q2;

q3 = Q + h .* (0.75 .* k2q);
p3 = P + h .* (0.75 .* k2p);
k3q = p3;
k3p = -omega^2 .* q3;

Q3 = Q + h .* ((2/9) .* k1q + (1/3) .* k2q + (4/9) .* k3q);
P3 = P + h .* ((2/9) .* k1p + (1/3) .* k2p + (4/9) .* k3p);
end

function [A_hist, t_hist, Q_hist, P_hist] = init_record(Q, P, n_steps, stride)
npts = numel(Q);
nrec = floor(n_steps / stride) + 1;
A_hist = zeros(nrec, 1);
t_hist = zeros(nrec, 1);
Q_hist = zeros(npts, nrec);
P_hist = zeros(npts, nrec);
A_hist(1) = cloud_area(Q, P);
Q_hist(:, 1) = Q(:);
P_hist(:, 1) = P(:);
end

function [A_hist, t_hist, Q_hist, P_hist] = save_record(A_hist, t_hist, Q_hist, P_hist, rec, t, Q, P)
A_hist(rec) = cloud_area(Q, P);
t_hist(rec) = t;
Q_hist(:, rec) = Q(:);
P_hist(:, rec) = P(:);
end

function A = cloud_area(Q, P)
k = convhull(Q(:), P(:));
A = polyarea(Q(k), P(k));
end
