function v_in = plane_wave_in(k, dx, dy, dz, x, y, z)
	v_in = exp(1i * k * (dx * x + dy * y + dz * z));
end
