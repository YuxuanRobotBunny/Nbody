function A = nbody_acceleration(R, params)
%NBODY_ACCELERATION Compute acceleration for each planet in an N-body system.
%
%   A = NBODY_ACCELERATION(R, params) returns an N-by-d array of accelerations.

params = nbody_prepare_params(params);

[N, d] = size(R);
if N ~= params.N || d ~= params.dim
    error('R must have size params.N-by-params.dim.');
end

G = params.G;
m0 = params.m0;
m = params.masses;
eps2 = params.eps_soft^2;

A = zeros(N, d);

% Force from the central star fixed at origin.
for i = 1:N
    ri = R(i, :);
    dist2 = sum(ri.^2) + eps2;
    if dist2 == 0
        error('Planet %d is at the origin. Central force is singular.', i);
    end
    dist3 = dist2 * sqrt(dist2);
    A(i, :) = A(i, :) - G * m0 * ri / dist3;
end

% Pairwise interaction between planets.
for i = 1:N-1
    for j = i+1:N
        rij = R(j, :) - R(i, :);
        dist2 = sum(rij.^2) + eps2;
        if dist2 == 0
            error('Planets %d and %d overlap. Pair force is singular.', i, j);
        end
        dist3 = dist2 * sqrt(dist2);
        a_ij = G * rij / dist3;
        A(i, :) = A(i, :) + m(j) * a_ij;
        A(j, :) = A(j, :) - m(i) * a_ij;
    end
end
end
