function everything = pro_open_everything(path, c, m, n)
	load(sprintf('%s/pro_%s.mat', path, generate_name(c, m, n)), '-mat');
end
