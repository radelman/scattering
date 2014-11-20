%
% Copyright (c) 2014, Ross Adelman, Nail A. Gumerov, and Ramani Duraiswami
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%

%
% Generate six figures.  The first three show a plane wave striking a prolate
% spheroid, oblate spheroid, and disk from an angle for k = 10.  The second
% three show the same, except for k = 25.
%
function generate_figure2abcdef()
	generate_pro_plane_wave_scat_figure(10.0, 1.0, pi - pi / 4.0, 1.5, 'hard', 5.0, 200);
	figure(3);
	export_fig(3, 'images/figure2a.pdf');
	close('all');
	generate_obl_plane_wave_scat_figure(10.0, 1.0, pi - pi / 4.0, 0.5, 'hard', 5.0, 200);
	figure(3);
	export_fig(3, 'images/figure2b.pdf');
	close('all');
	generate_obl_plane_wave_scat_figure(10.0, 1.0, pi - pi / 4.0, 0.0, 'hard', 5.0, 200);
	figure(3);
	export_fig(3, 'images/figure2c.pdf');
	close('all');
	generate_pro_plane_wave_scat_figure(25.0, 1.0, pi - pi / 4.0, 1.5, 'hard', 5.0, 200);
	figure(3);
	export_fig(3, 'images/figure2d.pdf');
	close('all');
	generate_obl_plane_wave_scat_figure(25.0, 1.0, pi - pi / 4.0, 0.5, 'hard', 5.0, 200);
	figure(3);
	export_fig(3, 'images/figure2e.pdf');
	close('all');
	generate_obl_plane_wave_scat_figure(25.0, 1.0, pi - pi / 4.0, 0.0, 'hard', 5.0, 200);
	figure(3);
	export_fig(3, 'images/figure2f.pdf');
	close('all');
end
