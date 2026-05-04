clear;
clc;

L = 60;
c = 1.0;
n_periods = 4;
N_list = [27, 81, 243];
methods = {'ETDRK2', 'SBDF2'};
dt_scan = [0.40, 0.30, 0.24, 0.18, 0.12, 0.09, 0.06, 0.045, 0.03, 0.0225];

fprintf('KdV empirical stability over %.1f periods\n', n_periods);
fprintf('method\t N\t dt_stable_max\t tested_dt_count\n');
for im = 1:numel(methods)
    method = methods{im};
    for iN = 1:numel(N_list)
        N = N_list(iN);
        grid = kdv_setup(N, L, struct('dealias', false));
        u0 = kdv_soliton(grid.x, 0.0, c, -15.0);
        u0_hat = fft(u0);
        T = n_periods * L / c;

        dt_stable = NaN;
        for j = 1:numel(dt_scan)
            dt = dt_scan(j);
            if run_stability_trial(method, grid, u0_hat, [0, T], dt)
                dt_stable = dt;
                break;
            end
        end
        fprintf('%s\t %d\t %.4f\t\t %d\n', method, N, dt_stable, numel(dt_scan));
    end
end

function ok = run_stability_trial(method, grid, u0_hat, tspan, dt)
try
    if strcmpi(method, 'ETDRK2')
        [~, uhat_hist, ~] = etdrk2_kdv(grid, tspan, u0_hat, dt);
    else
        [~, uhat_hist, ~] = sbdf2_kdv(grid, tspan, u0_hat, dt);
    end
    u = real(ifft(uhat_hist, [], 1));
    amp = max(abs(u(:)));
    ok = isfinite(amp) && amp < 5.0;
catch
    ok = false;
end
end
