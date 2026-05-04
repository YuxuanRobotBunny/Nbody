clear;
clc;

L = 60;
c = 1.0;
x0 = -15.0;
T = L / c;
N_list = [256, 128, 64, 32];
methods = {'ETDRK2', 'SBDF2'};

for im = 1:numel(methods)
    method = methods{im};
    figure('Name', sprintf('%s robustness', method));
    tiledlayout(2, 2, 'Padding', 'compact');
    fprintf('\nRobustness study for %s\n', method);
    fprintf('N\t dt_limit\t run_dt\t\t final_inf_error\n');

    for iN = 1:numel(N_list)
        N = N_list(iN);
        grid = kdv_setup(N, L, struct('dealias', false));
        u0 = kdv_soliton(grid.x, 0.0, c, x0);
        u0_hat = fft(u0);
        dt_limit = find_empirical_stability_limit(method, grid, u0_hat, T);
        run_dt = 0.5 * dt_limit;

        [~, uhat_hist] = run_method(method, grid, u0_hat, [0, T], run_dt);
        u_final = real(ifft(uhat_hist(:, end)));
        u_exact = kdv_soliton(grid.x, 0.0, c, x0);
        err_inf = norm(u_final - u_exact, inf);
        fprintf('%d\t %.4f\t\t %.4f\t %.3e\n', N, dt_limit, run_dt, err_inf);

        nexttile;
        plot(grid.x, u_exact, 'k--', 'LineWidth', 1.0, 'DisplayName', 'exact'); hold on;
        plot(grid.x, u_final, 'LineWidth', 1.2, 'DisplayName', 'numerical');
        grid on;
        xlabel('x');
        ylabel('u(x,T)');
        title(sprintf('N=%d, dt=%.4f', N, run_dt));
        legend('Location', 'best');
    end
end

function dt_limit = find_empirical_stability_limit(method, grid, u0_hat, T)
dt_scan = [0.40, 0.30, 0.24, 0.18, 0.12, 0.09, 0.06, 0.045, 0.03, 0.0225, 0.015];
dt_limit = dt_scan(end);
for j = 1:numel(dt_scan)
    dt = dt_scan(j);
    if run_stability_trial(method, grid, u0_hat, [0, T], dt)
        dt_limit = dt;
        return;
    end
end
end

function ok = run_stability_trial(method, grid, u0_hat, tspan, dt)
try
    [~, uhat_hist] = run_method(method, grid, u0_hat, tspan, dt);
    u = real(ifft(uhat_hist, [], 1));
    amp = max(abs(u(:)));
    ok = isfinite(amp) && amp < 5.0;
catch
    ok = false;
end
end

function [t, uhat_hist] = run_method(method, grid, u0_hat, tspan, dt)
if strcmpi(method, 'ETDRK2')
    [t, uhat_hist] = etdrk2_kdv(grid, tspan, u0_hat, dt);
else
    [t, uhat_hist] = sbdf2_kdv(grid, tspan, u0_hat, dt);
end
end
