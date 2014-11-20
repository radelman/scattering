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

function generate_data_for_all_figures()
	c = 10.0;
	for m = 0 : 29
		m
		
		for n = m : m + 29
			pro_calculate_everything(sprintf('%s/sphwv', pwd()), 2000, 500, c, m, n, 10, -300, 10, -300, 10, -300);
			everything = pro_load_everything('data', c, m, n, '');
			pro_save_everything('saved', everything);
			obl_calculate_everything(sprintf('%s/sphwv', pwd()), 2000, 500, c, m, n, 10, -300, 10, -300, 10, -300, 10, -300);
			everything = obl_load_everything('data', c, m, n, '');
			obl_save_everything('saved', everything);
		end
	end
	m = 10;
	for n = m : m + 29
		n
		
		pro_calculate_functions('', sprintf('%s/sphwv', pwd()), 2000, 300, c, m, n, 1.0 / 512.0, 'theta/pi', -7, 9.0, 1.0 / 64.0, 'x', 17);
		everything = pro_load_everything('data', c, m, n, '');
		pro_save_everything('saved', everything);
		obl_calculate_functions('', sprintf('%s/sphwv', pwd()), 2000, 300, c, m, n, 1.0 / 512.0, 'theta/pi', 8.0, 1.0 / 64.0, 'z', 17);
		everything = obl_load_everything('data', c, m, n, '');
		obl_save_everything('saved', everything);
	end
	c = 25.0;
	m = 49;
	for n = m : m + 49
		n
		
		obl_calculate_everything(sprintf('%s/sphwv', pwd()), 2000, 500, c, m, n, 10, -300, 10, -300, 10, -300, 10, -300);
		everything = obl_load_everything('data', c, m, n, '');
		obl_save_everything('saved', everything);
	end
end
