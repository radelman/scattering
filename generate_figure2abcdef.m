%
% Generate six figures.  The first three show a plane wave striking a prolate
% spheroid, oblate spheroid, and disk from an angle for k = 10.  The second
% three show the same, except for k = 25.
%
function generate_figure2abcdef()
	generate_pro_plane_wave_scat_figure(10.0, 1.0, pi - pi / 4.0, 1.5, 'hard', 5.0, 200);
	figure(3);
	set(gca, 'fontsize', 20);
	export_fig(3, 'images/figure2a.pdf');
	close('all');
	generate_obl_plane_wave_scat_figure(10.0, 1.0, pi - pi / 4.0, 0.5, 'hard', 5.0, 200);
	figure(3);
	set(gca, 'fontsize', 20);
	export_fig(3, 'images/figure2b.pdf');
	close('all');
	generate_obl_plane_wave_scat_figure(10.0, 1.0, pi - pi / 4.0, 0.0, 'hard', 5.0, 200);
	figure(3);
	set(gca, 'fontsize', 20);
	export_fig(3, 'images/figure2c.pdf');
	close('all');
	generate_pro_plane_wave_scat_figure(25.0, 1.0, pi - pi / 4.0, 1.5, 'hard', 5.0, 200);
	figure(3);
	set(gca, 'fontsize', 20);
	export_fig(3, 'images/figure2d.pdf');
	close('all');
	generate_obl_plane_wave_scat_figure(25.0, 1.0, pi - pi / 4.0, 0.5, 'hard', 5.0, 200);
	figure(3);
	set(gca, 'fontsize', 20);
	export_fig(3, 'images/figure2e.pdf');
	close('all');
	generate_obl_plane_wave_scat_figure(25.0, 1.0, pi - pi / 4.0, 0.0, 'hard', 5.0, 200);
	figure(3);
	set(gca, 'fontsize', 20);
	export_fig(3, 'images/figure2f.pdf');
	close('all');
end
