%
% Generate a test figure showing a plane wave striking an oblate spheroid from
% directly above for k = 10.
%
function generate_test_figure2()
	generate_obl_plane_wave_scat_figure(10.0, 1.0, pi, 0.5, 'hard', 5.0, 200);
	figure(3);
	set(gca, 'fontsize', 20);
	export_fig(3, 'images/test_figure2.pdf');
	close('all');
end
