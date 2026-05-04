clear;
clc;

% Single-planet circular orbit: fixed-step BS-RK23 verification.
params = struct('G', 1.0, 'm0', 1.0, 'masses', 1.0, 'dim', 2, 'eps_soft', 0.0);
params = nbody_prepare_params(params);

r0 = 1.0;
omega = 1.0;
T = 2 * pi / omega;

R0 = [r0, 0.0];
V0 = [0.0, omega * r0];
y0 = pack_state(R0, V0);

M_list = [40, 80, 160, 320, 640, 1280];
h_list = T ./ M_list;
err_max = zeros(size(M_list));

for k = 1:numel(M_list)
    h = h_list(k);
    [t, Y, ~] = bs23_fixed(@nbody_rhs, [0, T], y0, h, params);
    err_max(k) = max_position_error_vs_exact(t, Y, r0, omega, params);
end

p_est = log(err_max(1:end-1) ./ err_max(2:end)) ./ log(h_list(1:end-1) ./ h_list(2:end));
fprintf('Estimated order per refinement level:\n');
disp(p_est(:).');

figure;
loglog(h_list, err_max, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('h');
ylabel('max ||r_h(t)-r_{exact}(t)||_2');
title('BS-RK23 Fixed Step: Convergence for Circular Orbit');

target_err = 1e-2 * r0;
M_star = find_min_steps_for_error(y0, T, target_err, r0, omega, params);
h_star = T / M_star;
fprintf('Largest h with max position error <= %.3e is h = %.6e (M = %d)\n', target_err, h_star, M_star);

[t_c, Y_c, ~] = bs23_fixed(@nbody_rhs, [0, T], y0, h_star, params);
[t_f, Y_f, ~] = bs23_fixed(@nbody_rhs, [0, T], y0, h_star/2, params);

R_c = extract_planet_trajectory(Y_c, params, 1);
R_f = extract_planet_trajectory(Y_f, params, 1);
R_f_on_c = R_f(:, 1:2:end);
R_rich = R_f_on_c + (R_f_on_c - R_c) / (2^3 - 1);

[R_exact_c_x, R_exact_c_y] = exact_circle(t_c, r0, omega);
R_exact_c = [R_exact_c_x.'; R_exact_c_y.'];

err_c = vecnorm(R_c - R_exact_c, 2, 1);
err_f = vecnorm(R_f_on_c - R_exact_c, 2, 1);
err_r = vecnorm(R_rich - R_exact_c, 2, 1);

figure;
plot(t_c, err_c, 'LineWidth', 1.3); hold on;
plot(t_c, err_f, 'LineWidth', 1.3);
plot(t_c, err_r, 'LineWidth', 1.3);
grid on;
xlabel('t');
ylabel('position error');
legend('h', 'h/2 sampled on h grid', 'Richardson', 'Location', 'best');
title('Error vs Time (Linear Scale)');

figure;
semilogy(t_c, max(err_c, 1e-16), 'LineWidth', 1.3); hold on;
semilogy(t_c, max(err_f, 1e-16), 'LineWidth', 1.3);
semilogy(t_c, max(err_r, 1e-16), 'LineWidth', 1.3);
grid on;
xlabel('t');
ylabel('position error (log scale)');
legend('h', 'h/2 sampled on h grid', 'Richardson', 'Location', 'best');
title('Error vs Time (Log Scale)');

fprintf('Max errors: coarse = %.3e, fine = %.3e, richardson = %.3e\n', ...
    max(err_c), max(err_f), max(err_r));

function err_max = max_position_error_vs_exact(t, Y, r0, omega, params)
R = extract_planet_trajectory(Y, params, 1);
[x_ex, y_ex] = exact_circle(t, r0, omega);
R_ex = [x_ex.'; y_ex.'];
err = vecnorm(R - R_ex, 2, 1);
err_max = max(err);
end

function M_star = find_min_steps_for_error(y0, T, target_err, r0, omega, params)
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
    if M_hi > 2^20
        error('Failed to bracket target error for binary search.');
    end
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

M_star = M_hi;
end

function R = extract_planet_trajectory(Y, params, planet_idx)
R_all = reshape(Y(1:params.N*params.dim, :), params.N, params.dim, []);
R = squeeze(R_all(planet_idx, :, :));
end

function [x, y] = exact_circle(t, r0, omega)
x = r0 * cos(omega * t);
y = r0 * sin(omega * t);
end
