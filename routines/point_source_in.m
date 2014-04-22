function v_in = point_source_in(k, px, py, pz, x, y, z)
	r = sqrt((x - px) .^ 2 + (y - py) .^ 2 + (z - pz) .^ 2);
	v_in = exp(1i * k * r) ./ (4.0 * pi * r);
end
