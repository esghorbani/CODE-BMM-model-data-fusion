function phi = mode_shapesPP(x, n_modes, L)
%MODE_SHAPESPP Simply supported beam mode shapes (column vector)
%   x : scalar position
%   phi : (n_modes x 1)

if isempty(x) || ~isfinite(x)
    error('mode_shapesPP: x is empty or not finite.');
end

x = x(1);                 % enforce scalar
n = (1:n_modes).';        % column

phi = sqrt(2/L) * sin(n*pi*x/L);
end
