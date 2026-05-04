clear;
clc;

L = 60;
c = 1.0;
x0 = -15.0;
N = 256;
grid = kdv_setup(N, L, struct('dealias', false));
u0 = kdv_soliton(grid.x, 0.0, c, x0);
u0_hat = fft(u0);

T_full = L / c;
T = T_full / 4;
dt_list = [0.12, 0.06, 0.03, 0.015];
methods = {'ETDRK2', 'SBDF2'};

fprintf('KdV accuracy study on T = %.3f (one quarter period)\n', T);
for im = 1:numel(methods)
    method = methods{im};
    fprintf('\nMethod: %s\n', method);
    err_final_refine = zeros(numel(dt_list) - 1, 1);
    err_final_exact = zeros(numel(dt_list), 1);

    cache = cell(numel(dt_list), 2);
    for j = 1:numel(dt_list)
        dt = dt_list(j);
        [t, uhat_hist] = run_method(method, grid, u0_hat, [0, T], dt);
        u_hist = real(ifft(uhat_hist, [], 1));
        u_exact = kdv_soliton(grid.x, t(end), c, x0);
        err_final_exact(j) = norm(u_hist(:, end) - u_exact, inf);
        cache{j, 1} = t;
        cache{j, 2} = u_hist;
    end

    for j = 1:(numel(dt_list) - 1)
        u_coarse = cache{j, 2};
        u_fine = cache{j + 1, 2};
        u_fine_on_coarse = u_fine(:, 1:2:end);
        err_series = vecnorm(u_coarse - u_fine_on_coarse, inf, 1);
        err_final_refine(j) = err_series(end);

        figure('Name', sprintf('%s refine dt=%.4f', method, dt_list(j)));
        semilogy(cache{j, 1}, max(err_series, 1e-16), 'LineWidth', 1.3);
        grid on;
        xlabel('t');
        ylabel('||u_{\Delta t}-u_{\Delta t/2}||_{\infty}');
        title(sprintf('%s successive refinement, dt=%.4f', method, dt_list(j)));
    end

    p_refine = log(err_final_refine(1:end-1) ./ err_final_refine(2:end)) ./ log(2);
    p_exact = log(err_final_exact(1:end-1) ./ err_final_exact(2:end)) ./ log(2);
    disp(table(dt_list(1:end-1).', err_final_refine, [p_refine; NaN], ...
        'VariableNames', {'dt', 'successive_refinement_error', 'order_est'}));
    disp(table(dt_list.', err_final_exact, [p_exact; NaN], ...
        'VariableNames', {'dt', 'exact_error_inf', 'order_est'}));

    figure('Name', sprintf('%s exact final error', method));
    loglog(dt_list, err_final_exact, 'o-', 'LineWidth', 1.4);
    grid on;
    xlabel('\Deltat');
    ylabel('||u(\cdot,T)-u_{exact}(\cdot,T)||_{\infty}');
    title(sprintf('%s PDE accuracy at final time', method));
end

function [t, uhat_hist] = run_method(method, grid, u0_hat, tspan, dt)
if strcmpi(method, 'ETDRK2')
    [t, uhat_hist] = etdrk2_kdv(grid, tspan, u0_hat, dt);
else
    [t, uhat_hist] = sbdf2_kdv(grid, tspan, u0_hat, dt);
end
end
