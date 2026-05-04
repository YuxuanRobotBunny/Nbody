function params = nbody_prepare_params(params)
%NBODY_PREPARE_PARAMS Validate and fill defaults for N-body parameters.

if ~isfield(params, 'masses')
    error('params.masses is required.');
end
params.masses = params.masses(:);

if ~isfield(params, 'N')
    params.N = numel(params.masses);
end
if params.N ~= numel(params.masses)
    error('params.N must match numel(params.masses).');
end

if ~isfield(params, 'dim')
    params.dim = 2;
end
if ~(params.dim == 2 || params.dim == 3)
    error('params.dim must be 2 or 3.');
end

if ~isfield(params, 'G')
    params.G = 1.0;
end
if ~isfield(params, 'm0')
    params.m0 = 1.0;
end
if ~isfield(params, 'eps_soft')
    params.eps_soft = 0.0;
end
if params.eps_soft < 0
    error('params.eps_soft must be nonnegative.');
end
end
