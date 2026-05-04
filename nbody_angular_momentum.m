function L = nbody_angular_momentum(y, params)
%NBODY_ANGULAR_MOMENTUM Angular momentum of the planets.
%
%   In 2D this returns scalar Lz.
%   In 3D this returns a 1-by-3 vector.

params = nbody_prepare_params(params);
[R, V] = unpack_state(y, params.N, params.dim);
m = params.masses;

if params.dim == 2
    L = sum(m .* (R(:, 1) .* V(:, 2) - R(:, 2) .* V(:, 1)));
    return;
end

L = zeros(1, 3);
for i = 1:params.N
    L = L + m(i) * cross(R(i, :), V(i, :));
end
end
