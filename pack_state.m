function y = pack_state(R, V)
%PACK_STATE Flatten position and velocity arrays into a single state vector.
%
%   y = PACK_STATE(R, V) returns y = [R(:); V(:)] where R and V are N-by-d.

if ~isequal(size(R), size(V))
    error('R and V must have identical sizes.');
end

y = [R(:); V(:)];
end
