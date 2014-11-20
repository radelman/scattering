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
% Calculate and save the oblate spheroidal wave functions.  This function takes
% care of calling obl_sphwv for you.  It also hides some of the features of
% obl_sphwv, so there is a trade-off.  Because the output from obl_sphwv is in
% ASCII and, thus, takes up a lot of disk space, this function ZIPs up
% everything at the end.
%
% Arguments:
%     name - a name for the set of oblate spheroidal wave functions about to
%            be calculated.  at the end, when everything is ZIP'd up, this name
%            is appended to the end of that ZIP.  can be left blank ('')
%     path - the directory in which obl_sphwv is located
%     max_memory - the maximum amount of memory, in MB, that obl_sphwv can use
%                  before automatically terminating
%     precision - the number of bits of precision obl_sphwv should use
%     c - k * a, where k is the wavenumber and 2a is the interfocal distance
%     m - one modal value
%     n - the other modal value, which is equal to m, m + 1, ...
%     d_S1 - the spacing between argument values when calculating the oblate
%            spheroidal angle functions.
%     arg_type_S1 - whether the argument for the oblate spheroidal angle
%                   functions should be eta (-1 to 1) or theta (0 to pi).  in
%                   the latter, the argument really goes from 0 to 1, and
%                   internally, obl_sphwv multiplies this by pi.  this is
%                   because working with fractions of pi can be messy.
%                   depending on whether you enter 'eta' or 'theta/pi', you
%                   should pick your choice of d_S1 accordingly
%     b - the range of values over which to calculate the oblate spheroidal
%         radial functions.  the range will be [0, b]
%     d_R - the spacing between argument values when calculating the oblate
%           spheroidal radial functions.
%     arg_type_R - whether the argument for the oblate spheroidal radial
%                  functions should be xi (>= 0) or z (>= 0).  in the latter,
%                  obl_sphwv internally transforms z into xi using xi = z.
%                  in other words, they're equivalent.  this option is here so
%                  that this function matches up to its prolate counterpart.
%     p - how many digits of precision obl_sphwv should output
% Return Values:
%     None.
%
function obl_calculate_functions(name, path, max_memory, precision, c, m, n, d_S1, arg_type_S1, b, d_R, arg_type_R, p)
	unzip(sprintf('data/obl_%s.zip', generate_name(c, m, n)), 'data');
	fprintf('calculating S1...\n');
	if (strcmp(arg_type_S1, 'eta'))
		command = sprintf('"%s/obl_sphwv" -max_memory %d -precision %d -verbose n -c %s -m %d -n %d -w S1 -a -1.0 -b 1.0 -d %s -arg_type eta -p %d > data/obl_%s_S1.txt', path, max_memory, precision, nice_number(c, 20), m, n, nice_number(d_S1, 20), p, generate_name(c, m, n));
	else
		command = sprintf('"%s/obl_sphwv" -max_memory %d -precision %d -verbose n -c %s -m %d -n %d -w S1 -a 0.0 -b 1.0 -d %s -arg_type theta/pi -p %d > data/obl_%s_S1.txt', path, max_memory, precision, nice_number(c, 20), m, n, nice_number(d_S1, 20), p, generate_name(c, m, n));
	end
	fprintf('%s\n', command);
	[ ...
	status, result ...
	] = system(command, '-echo');
	if (status == 1)
		return;
	end
	fprintf('calculating R...\n');
	if (strcmp(arg_type_R, 'xi'))
		command = sprintf('"%s/obl_sphwv" -max_memory %d -precision %d -verbose n -c %s -m %d -n %d -w R -a 0.0 -b %s -d %s -arg_type xi -which R1_1,R1_2,R2_1,R2_2,R2_31,R2_32 -p %d > data/obl_%s_R.txt', path, max_memory, precision, nice_number(c, 20), m, n, nice_number(b, 20), nice_number(d_R, 20), p, generate_name(c, m, n));
	else
		command = sprintf('"%s/obl_sphwv" -max_memory %d -precision %d -verbose n -c %s -m %d -n %d -w R -a 0.0 -b %s -d %s -arg_type z -which R1_1,R1_2,R2_1,R2_2,R2_31,R2_32 -p %d > data/obl_%s_R.txt', path, max_memory, precision, nice_number(c, 20), m, n, nice_number(b, 20), nice_number(d_R, 20), p, generate_name(c, m, n));
	end
	fprintf('%s\n', command);
	[ ...
	status, result ...
	] = system(command, '-echo');
	if (status == 1)
		return;
	end
	files = {sprintf('obl_%s_S1.txt', generate_name(c, m, n)), sprintf('obl_%s_R.txt', generate_name(c, m, n))};
	zip(sprintf('data/obl_%s_%s.zip', generate_name(c, m, n), name), files, 'data');
	delete(sprintf('data/obl_%s_*.txt', generate_name(c, m, n)));
end
