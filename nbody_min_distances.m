function [d_star_min, d_pair_min] = nbody_min_distances(y, params)
%NBODY_MIN_DISTANCES Minimum distance to star and minimum pair distance.

params = nbody_prepare_params(params);
[R, ~] = unpack_state(y, params.N, params.dim);

dist_star = sqrt(sum(R.^2, 2));
d_star_min = min(dist_star);

d_pair_min = inf;
for i = 1:params.N-1
    for j = i+1:params.N
        dij = norm(R(j, :) - R(i, :), 2);
        if dij < d_pair_min
            d_pair_min = dij;
        end
    end
end
end
