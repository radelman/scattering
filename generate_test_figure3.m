%
% Generate a test figure showing a plane wave striking a disk from directly
% above for k = 10.
%
function generate_test_figure3()
	generate_obl_plane_wave_scat_figure(10.0, 1.0, pi, 0.0, 'hard', 5.0, 200);
	figure(3);
	set(gca, 'fontsize', 20);
	export_fig(3, 'images/test_figure3.pdf');
	close('all');
end
