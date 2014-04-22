generate_figure1();
export_fig(1, 'images/figure1.pdf');
close('all');

generate_figure2(10.0, 'Sound-Hard Prolate Spheroid, k = 10');
export_fig(1, 'images/figure2a.pdf');
close('all');
generate_figure2(25.0, 'Sound-Hard Prolate Spheroid, k = 25');
export_fig(1, 'images/figure2b.pdf');
close('all');

generate_figure3();
export_fig(1, 'images/figure3.pdf');
close('all');

generate_figure4(10.0, 0.5, 'Sound-Hard Oblate Spheroid, k = 10');
export_fig(1, 'images/figure4a.pdf');
close('all');
generate_figure4(25.0, 0.5, 'Sound-Hard Oblate Spheroid, k = 25');
export_fig(1, 'images/figure4b.pdf');
close('all');
generate_figure4(10.0, 0.0, 'Sound-Hard Disk, k = 10');
export_fig(1, 'images/figure4c.pdf');
close('all');
generate_figure4(25.0, 0.0, 'Sound-Hard Disk, k = 25');
export_fig(1, 'images/figure4d.pdf');
close('all');

generate_figure6(10.0, 'Sound-Hard Disk, k = 10');
export_fig(1, 'images/figure6a.pdf');
close('all');
generate_figure6(25.0, 'Sound-Hard Disk, k = 25');
export_fig(1, 'images/figure6b.pdf');
close('all');

generate_figure7(10.0, 'Sound-Hard Disk, k = 10');
export_fig(1, 'images/figure7a.pdf');
close('all');
generate_figure7(25.0, 'Sound-Hard Disk, k = 25');
export_fig(1, 'images/figure7b.pdf');
close('all');
generate_figure7(100.0, 'Sound-Hard Disk, k = 100');
export_fig(1, 'images/figure7c.pdf');
close('all');
generate_figure7(250.0, 'Sound-Hard Disk, k = 250');
export_fig(1, 'images/figure7d.pdf');
close('all');
