generate_test_figure1(10.0, 'Sound-Hard Prolate Spheroid, k = 10');
export_fig(1, 'images/test_figure1a.pdf');
close('all');
generate_test_figure1(25.0, 'Sound-Hard Prolate Spheroid, k = 25');
export_fig(1, 'images/test_figure1b.pdf');
close('all');

generate_test_figure2(10.0, 0.5, 'Sound-Hard Oblate Spheroid, k = 10');
export_fig(1, 'images/test_figure2a.pdf');
close('all');
generate_test_figure2(25.0, 0.5, 'Sound-Hard Oblate Spheroid, k = 25');
export_fig(1, 'images/test_figure2b.pdf');
close('all');
generate_test_figure2(10.0, 0.0, 'Sound-Hard Disk, k = 10');
export_fig(1, 'images/test_figure2c.pdf');
close('all');
generate_test_figure2(25.0, 0.0, 'Sound-Hard Disk, k = 25');
export_fig(1, 'images/test_figure2d.pdf');
close('all');
