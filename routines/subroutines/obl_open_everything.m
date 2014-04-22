function everything = obl_open_everything(path, c, m, n)
	load(sprintf('%s/obl_%s.mat', path, generate_name(c, m, n)), '-mat');
end
