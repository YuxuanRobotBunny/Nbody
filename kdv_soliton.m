function u = kdv_soliton(x, t, c, x0)
%KDV_SOLITON Traveling-wave soliton profile for KdV on the line.

if nargin < 4
    x0 = 0.0;
end
xi = x - x0 - c * t;
u = 0.5 * c * sech(0.5 * sqrt(c) * xi).^2;
end
