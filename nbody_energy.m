function E = nbody_energy(y, params)
%NBODY_ENERGY Total mechanical energy of the N-body system.
%
%   E = NBODY_ENERGY(y, params) returns kinetic + potential energy.

params = nbody_prepare_params(params);
[R, V] = unpack_state(y, params.N, params.dim);

G = params.G;
m0 = params.m0;
m = params.masses;
eps2 = params.eps_soft^2;

kin = 0.5 * sum(m .* sum(V.^2, 2));

% Potential due to central star.
r_norm = sqrt(sum(R.^2, 2) + eps2);
pot_star = -G * m0 * sum(m ./ r_norm);

% Pairwise potential.
pot_pair = 0.0;
for i = 1:params.N-1
    for j = i+1:params.N
        rij = R(j, :) - R(i, :);
        dist = sqrt(sum(rij.^2) + eps2);
        pot_pair = pot_pair - G * m(i) * m(j) / dist;
    end
end

E = kin + pot_star + pot_pair;
end
