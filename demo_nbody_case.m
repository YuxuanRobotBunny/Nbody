clear;
clc;

% N-body example for part 2.2 with diagnostics and optional movie output.
params = struct();
params.G = 1.0;
params.m0 = 1.0;
params.masses = [3e-4; 4e-4; 5e-4; 6e-4; 7e-4];
params.dim = 2;
params.eps_soft = 0.0;
params = nbody_prepare_params(params);

N = params.N;
mu = params.G * params.m0;

radii = [1.00; 1.35; 1.75; 2.20; 2.70];
angles = [0.10; 1.15; 2.05; 3.10; 4.20];
speed_scale = [1.00; 0.97; 1.04; 0.94; 1.06];

R0 = [radii .* cos(angles), radii .* sin(angles)];
V0 = zeros(N, 2);
for i = 1:N
    tangent = [-sin(angles(i)), cos(angles(i))];
    v_circ = sqrt(mu / radii(i));
    V0(i, :) = speed_scale(i) * v_circ * tangent;
end

% Mild perturbation to increase interactions without immediate close-collision.
V0(2, :) = V0(2, :) + [0.00, -0.02];
V0(4, :) = V0(4, :) + [0.01, 0.00];

y0 = pack_state(R0, V0);
t_end = 350.0;
h0 = 0.02;
tol = 1e-6;

opts = struct( ...
    'err_mode', 'full', ...
    'rtol', tol, ...
    'atol', 1e-10, ...
    'controller', 'pi', ...
    'fac_max', 2.0, ...
    'fac_max_reject', 0.8, ...
    'h_min', 1e-7, ...
    'h_max', 0.2);
[t, Y, stats] = bs23_adaptive(@nbody_rhs, [0, t_end], y0, h0, tol, params, opts);

fprintf('Adaptive run complete: accepted=%d, rejected=%d, nfev=%d\n', ...
    stats.accepted_steps, stats.rejected_steps, stats.nfev);

R = reshape(Y(1:N*params.dim, :), N, params.dim, []);

figure;
hold on;
plot(0, 0, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 10, 'DisplayName', 'star');
for i = 1:N
    traj = squeeze(R(i, :, :));
    plot(traj(1, :), traj(2, :), 'LineWidth', 1.2, 'DisplayName', sprintf('planet %d', i));
end
axis equal;
grid on;
xlabel('x');
ylabel('y');
title('N-body trajectories (adaptive BS-RK23)');
legend('Location', 'bestoutside');

E = zeros(numel(t), 1);
L = zeros(numel(t), 1);
d_star_min = inf;
d_pair_min = inf;
for k = 1:numel(t)
    E(k) = nbody_energy(Y(:, k), params);
    L(k) = nbody_angular_momentum(Y(:, k), params);
    [d_star_k, d_pair_k] = nbody_min_distances(Y(:, k), params);
    d_star_min = min(d_star_min, d_star_k);
    d_pair_min = min(d_pair_min, d_pair_k);
end

fprintf('Minimum distance to star: %.6e\n', d_star_min);
fprintf('Minimum inter-planet distance: %.6e\n', d_pair_min);

dE = E - E(1);
dL = L - L(1);
fprintf('Energy drift range: [%.6e, %.6e], span=%.6e\n', min(dE), max(dE), max(dE)-min(dE));
fprintf('Angular momentum drift range: [%.6e, %.6e], span=%.6e\n', ...
    min(dL), max(dL), max(dL)-min(dL));

figure;
plot(t, dE, 'LineWidth', 1.2);
grid on;
xlabel('t');
ylabel('\Delta E(t) = E(t)-E(0)');
title('Energy drift diagnostic');

figure;
plot(t, dL, 'LineWidth', 1.2);
grid on;
xlabel('t');
ylabel('\Delta L_z(t) = L_z(t)-L_z(0)');
title('Angular momentum drift diagnostic');

figure;
acc = stats.accepted_attempt;
rej = ~acc;
plot(stats.t_attempt(acc), stats.h_attempt(acc), 'b.-', 'DisplayName', 'accepted'); hold on;
plot(stats.t_attempt(rej), stats.h_attempt(rej), 'rx', 'DisplayName', 'rejected');
grid on;
xlabel('t');
ylabel('attempted h');
title('Adaptive step size history');
legend('Location', 'best');

write_movie = false;  % Set true if your MATLAB installation supports MPEG-4 profile.
movie_name = 'nbody_demo.mp4';
if write_movie
    try
        make_movie(t, R, movie_name);
        fprintf('Movie saved to %s\n', movie_name);
    catch ME
        warning('Movie creation failed: %s', ME.message);
    end
end

function make_movie(t, R, movie_name)
[N, ~, nt] = size(R);
x_all = reshape(R(:, 1, :), N, nt);
y_all = reshape(R(:, 2, :), N, nt);

lim = 1.1 * max(max(abs([x_all(:), y_all(:)])));

v = VideoWriter(movie_name, 'MPEG-4');
v.FrameRate = 30;
open(v);

fig = figure('Color', 'w');
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
xlim(ax, [-lim, lim]);
ylim(ax, [-lim, lim]);
xlabel(ax, 'x');
ylabel(ax, 'y');
title(ax, 'N-body dynamics');
plot(ax, 0, 0, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 11);

colors = lines(N);
trail = gobjects(N, 1);
body = gobjects(N, 1);
for i = 1:N
    trail(i) = plot(ax, x_all(i, 1), y_all(i, 1), '-', 'Color', colors(i, :), 'LineWidth', 1.0);
    body(i) = plot(ax, x_all(i, 1), y_all(i, 1), 'o', 'Color', colors(i, :), ...
        'MarkerFaceColor', colors(i, :), 'MarkerSize', 5);
end
time_text = text(ax, -0.95*lim, 0.92*lim, sprintf('t = %.2f', t(1)));

stride = max(1, floor(nt / 1500));
for k = 1:stride:nt
    for i = 1:N
        set(trail(i), 'XData', x_all(i, 1:k), 'YData', y_all(i, 1:k));
        set(body(i), 'XData', x_all(i, k), 'YData', y_all(i, k));
    end
    set(time_text, 'String', sprintf('t = %.2f', t(k)));
    drawnow;
    writeVideo(v, getframe(fig));
end

close(v);
close(fig);
end
