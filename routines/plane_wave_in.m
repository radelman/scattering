%
% Calculate the incident field due to a plane wave.
%
% Arguments:
%     k - the wavenumber
%     dx, dy, dz - the direction of the plane wave in Cartesian coordinates
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
% Return Values:
%     v_in - the incident field
%     grad_in_cart - the gradient of the incident field in Cartesian
%                    coordinates
%
function [v_in, grad_in_cart] = plane_wave_in(k, dx, dy, dz, x, y, z)
	v_in = exp(1i * k * (dx * x + dy * y + dz * z));
	grad_in_cart = 1i * k * [dx; dy; dz] * v_in;
end
