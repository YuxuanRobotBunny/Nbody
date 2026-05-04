clear;
clc;

L = 80;
N = 512;
grid = kdv_setup(N, L, struct('dealias', false));

c1 = 1.0;
c2 = 2.0;
x1 = -20.0;
x2 = 8.0;
u0 = kdv_soliton(grid.x, 0.0, c1, x1) + kdv_soliton(grid.x, 0.0, c2, x2);
u0_hat = fft(u0);

T = 22.0;
dt = 0.015;
[t, uhat_hist, stats] = etdrk2_kdv(grid, [0, T], u0_hat, dt);
u_hist = real(ifft(uhat_hist, [], 1));

fprintf('Two-soliton ETDRK2 run: N=%d, dt=%.4f, nonlinear evals=%d\n', ...
    N, dt, stats.nonlinear_evals);

figure('Name', 'Two-soliton snapshots');
snap_idx = round(linspace(1, numel(t), 4));
for j = 1:numel(snap_idx)
    subplot(2, 2, j);
    idx = snap_idx(j);
    plot(grid.x, u_hist(:, idx), 'LineWidth', 1.3);
    grid on;
    xlabel('x');
    ylabel('u');
    title(sprintf('t = %.2f', t(idx)));
end

write_movie = false;
movie_name = 'kdv_two_solitons.mp4';
if write_movie
    make_movie(grid.x, t, u_hist, movie_name);
end

function make_movie(x, t, u_hist, movie_name)
v = VideoWriter(movie_name, 'MPEG-4');
v.FrameRate = 24;
open(v);

fig = figure('Color', 'w');
ax = axes(fig);
ylim(ax, [min(u_hist(:)) - 0.05, max(u_hist(:)) + 0.05]);
xlim(ax, [min(x), max(x)]);
grid(ax, 'on');
xlabel(ax, 'x');
ylabel(ax, 'u');
title(ax, 'KdV two-soliton interaction');
curve = plot(ax, x, u_hist(:, 1), 'LineWidth', 1.5);
time_text = text(ax, x(5), max(u_hist(:)), sprintf('t = %.2f', t(1)));

stride = max(1, floor(numel(t) / 1200));
for k = 1:stride:numel(t)
    set(curve, 'YData', u_hist(:, k));
    set(time_text, 'String', sprintf('t = %.2f', t(k)));
    drawnow;
    writeVideo(v, getframe(fig));
end

close(v);
close(fig);
end
