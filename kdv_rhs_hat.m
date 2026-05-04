function dudt_hat = kdv_rhs_hat(uhat, grid)
%KDV_RHS_HAT Right-hand side for Fourier-space KdV ODEs.

dudt_hat = grid.lin .* uhat + kdv_nonlinear_hat(uhat, grid);
dudt_hat(grid.mean_idx) = 0.0;
if grid.has_unmatched_mode
    dudt_hat(grid.unmatched_idx) = 0.0;
end
end
