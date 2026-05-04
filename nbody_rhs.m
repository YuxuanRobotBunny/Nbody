function dydt = nbody_rhs(~, y, params)
%NBODY_RHS First-order form of the N-body equations.
%
%   dydt = NBODY_RHS(t, y, params) returns [V(:); A(:)] where A is computed
%   from Newtonian gravity with a fixed central star at the origin.

params = nbody_prepare_params(params);
[R, V] = unpack_state(y, params.N, params.dim);
A = nbody_acceleration(R, params);
dydt = [V(:); A(:)];
end
