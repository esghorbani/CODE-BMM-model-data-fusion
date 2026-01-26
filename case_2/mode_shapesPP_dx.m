function phi_x = mode_shapesPP_dx(x, n_modes, L)
%MODE_SHAPESPP_DX First derivative of simply supported beam mode shapes
%   x : scalar position
%   phi_x : (n_modes x 1)

if isempty(x) || ~isfinite(x)
    error('mode_shapesPP_dx: x is empty or not finite.');
end

x = x(1);                 % enforce scalar
n = (1:n_modes).';        % column

phi_x = sqrt(2/L) * (n*pi/L) .* cos(n*pi*x/L);
end
