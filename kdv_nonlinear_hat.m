function Bhat = kdv_nonlinear_hat(uhat, grid)
%KDV_NONLINEAR_HAT Fourier coefficients of -3 * d_x (u^2).

u = real(ifft(uhat));
u2_hat = fft(u.^2);
if grid.dealias
    u2_hat = dealias_hat(u2_hat, grid);
end
Bhat = -3i * grid.k .* u2_hat;
Bhat(grid.mean_idx) = 0.0;
if grid.has_unmatched_mode
    Bhat(grid.unmatched_idx) = 0.0;
end
end

function uhat = dealias_hat(uhat, grid)
N = grid.N;
kc = floor(N / 3);
mask = false(N, 1);
if mod(N, 2) == 0
    mask(1:(kc + 1)) = true;
    mask((N - kc + 1):N) = true;
else
    mask(1:(kc + 1)) = true;
    mask((N - kc + 1):N) = true;
end
uhat(~mask) = 0.0;
end
