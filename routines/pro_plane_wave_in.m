%
% Calculate the incident field due to a plane wave using prolate spheroidal
% wave functions.
%
% Arguments:
%     k - the wavenumber
%     a - the interfocal distance divided by two
%     theta0 - the direction of the plane wave in spherical coordinates
%     path - the directory in which the prolate spheroidal wave function have
%            been precomputed and stored
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
% Return Values:
%     v_in - the incident field
%     max_abs_change - the maximum absolute change per mode over all the
%                      evaluation points
%
function [v_in, max_abs_change] = pro_plane_wave_in(k, a, theta0, path, x, y, z)
	[ ...
	v_in, max_abs_change ...
	] = pro_calculate_sum(k, a, path, x, y, z, @(k, a, c, m, n, everything, eta, xi, phi)(calculate_plane_wave_in_term(k, a, c, m, n, everything, eta, xi, phi, theta0)));
end
