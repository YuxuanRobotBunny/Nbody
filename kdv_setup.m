function grid = kdv_setup(N, L, opts)
%KDV_SETUP Build Fourier grid data for periodic KdV on [-L/2, L/2].

if nargin < 3 || isempty(opts)
    opts = struct();
end
if mod(N, 1) ~= 0 || N < 8
    error('N must be an integer >= 8.');
end
if L <= 0
    error('L must be positive.');
end

dealias = get_opt(opts, 'dealias', false);

x = linspace(-L/2, L/2, N + 1).';
x(end) = [];

if mod(N, 2) == 0
    k = (2 * pi / L) * [0:(N/2 - 1), 0, (-N/2 + 1):-1].';
else
    k = (2 * pi / L) * [0:((N - 1) / 2), -((N - 1) / 2):-1].';
end

grid = struct();
grid.N = N;
grid.L = L;
grid.x = x;
grid.k = k;
grid.lin = 1i * k.^3;
grid.dealias = dealias;
grid.mean_idx = 1;
grid.has_unmatched_mode = mod(N, 2) == 0;
if grid.has_unmatched_mode
    grid.unmatched_idx = N / 2 + 1;
else
    grid.unmatched_idx = [];
end
end

function v = get_opt(opts, name, default_value)
if isfield(opts, name)
    v = opts.(name);
else
    v = default_value;
end
end
