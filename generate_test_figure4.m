%
% Generate a test figure showing an incident field due to a point source
% striking a disk from directly above for k = 10.
%
function generate_test_figure4()
	generate_obl_point_source_scat_figure(10.0, 1.0, 0.0, 3.0, 0.0, 'hard', 5.0, 200);
	figure(3);
	set(gca, 'fontsize', 20);
	export_fig(3, 'images/test_figure4.pdf');
	close('all');
end
