function [R, V] = unpack_state(y, N, d)
%UNPACK_STATE Recover position and velocity arrays from a state vector.
%
%   [R, V] = UNPACK_STATE(y, N, d) reshapes y into R and V with size N-by-d.

expected_len = 2 * N * d;
if numel(y) ~= expected_len
    error('State vector length mismatch: expected %d, got %d.', expected_len, numel(y));
end

R = reshape(y(1:N*d), N, d);
V = reshape(y(N*d+1:end), N, d);
end
