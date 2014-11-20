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
% Calculate and save the prolate spheroidal wave functions.  This function
% takes care of calling pro_sphwv for you.  It also hides some of the
% features of pro_sphwv, so there is a trade-off.  Because the output from
% pro_sphwv is in ASCII and, thus, takes up a lot of disk space, this function
% ZIPs up everything at the end.
%
% Arguments:
%     name - a name for the set of prolate spheroidal wave functions about to
%            be calculated.  at the end, when everything is ZIP'd up, this name
%            is appended to the end of that ZIP.  can be left blank ('')
%     path - the directory in which pro_sphwv is located
%     max_memory - the maximum amount of memory, in MB, that pro_sphwv can use
%                  before automatically terminating
%     precision - the number of bits of precision pro_sphwv should use
%     c - k * a, where k is the wavenumber and 2a is the interfocal distance
%     m - one modal value
%     n - the other modal value, which is equal to m, m + 1, ...
%     d_S1 - the spacing between argument values when calculating the prolate
%            spheroidal angle functions.
%     arg_type_S1 - whether the argument for the prolate spheroidal angle
%                   functions should be eta (-1 to 1) or theta (0 to pi).  in
%                   the latter, the argument really goes from 0 to 1, and
%                   internally, pro_sphwv multiplies this by pi.  this is
%                   because working with fractions of pi can be messy.
%                   depending on whether you enter 'eta' or 'theta/pi', you
%                   should pick your choice of d_S1 accordingly
%     e_min - in addition to calculating the prolate spheroidal radial
%             functions from 1 to b in increments of d_R, pro_sphwv also
%             calculates the values, 2 ^ e_min < 2 ^ (e_min + 1) < ..., up to,
%             but not including, d_R.  make sure to choose e_min so that
%             there's at least one of these values being calculated
%     b - the range of values over which to calculate the prolate spheroidal
%         radial functions.  the range will be [1, b]
%     d_R - the spacing between argument values when calculating the prolate
%           spheroidal radial functions.
%     arg_type_R - whether the argument for the prolate spheroidal radial
%                  functions should be xi (>= 1) or x (>= 0).  in the latter,
%                  pro_sphwv internally transforms x into xi using
%                  sqrt(x ^ 2 + 1).  depending on whether you enter 'xi' or
%                  'x', you should pick your choice of d_R accordingly
%     p - how many digits of precision pro_sphwv should output
% Return Values:
%     None.
%
function pro_calculate_functions(name, path, max_memory, precision, c, m, n, d_S1, arg_type_S1, e_min, b, d_R, arg_type_R, p)
	unzip(sprintf('data/pro_%s.zip', generate_name(c, m, n)), 'data');
	fprintf('calculating S1...\n');
	if (strcmp(arg_type_S1, 'eta'))
		command = sprintf('"%s/pro_sphwv" -max_memory %d -precision %d -verbose n -c %s -m %d -n %d -w S1 -a -1.0 -b 1.0 -d %s -arg_type eta -p %d > data/pro_%s_S1.txt', path, max_memory, precision, nice_number(c, 20), m, n, nice_number(d_S1, 20), p, generate_name(c, m, n));
	else
		command = sprintf('"%s/pro_sphwv" -max_memory %d -precision %d -verbose n -c %s -m %d -n %d -w S1 -a 0.0 -b 1.0 -d %s -arg_type theta/pi -p %d > data/pro_%s_S1.txt', path, max_memory, precision, nice_number(c, 20), m, n, nice_number(d_S1, 20), p, generate_name(c, m, n));
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
		command = sprintf('"%s/pro_sphwv" -max_memory %d -precision %d -verbose n -c %s -m %d -n %d -w R -a 1.0 -b %s -d %s -arg_type xi -which R1_1,R1_2,R2_1,R2_2 -p %d > data/pro_%s_R.txt', path, max_memory, precision, nice_number(c, 20), m, n, nice_number(b, 20), nice_number(d_R, 20), p, generate_name(c, m, n));
	else
		command = sprintf('"%s/pro_sphwv" -max_memory %d -precision %d -verbose n -c %s -m %d -n %d -w R -a 0.0 -b %s -d %s -arg_type x -which R1_1,R1_2,R2_1,R2_2 -p %d > data/pro_%s_R.txt', path, max_memory, precision, nice_number(c, 20), m, n, nice_number(b, 20), nice_number(d_R, 20), p, generate_name(c, m, n));
	end
	fprintf('%s\n', command);
	[ ...
	status, result ...
	] = system(command, '-echo');
	if (status == 1)
		return;
	end
	for e = e_min : log2(d_R) - 1
		if (e > e_min)
			output_type = '>>';
		else
			output_type = '>';
		end
		if (strcmp(arg_type_R, 'xi'))
			command = sprintf('"%s/pro_sphwv" -max_memory %d -precision %d -verbose n -c %s -m %d -n %d -w R -a "1+2^%d" -b "1+2^%d" -d 1.0 -arg_type xi -which R1_1,R1_2,R2_1,R2_2 -p %d %s data/pro_%s_R_small.txt', path, max_memory, precision, nice_number(c, 20), m, n, e, e, p, output_type, generate_name(c, m, n));
		else
			command = sprintf('"%s/pro_sphwv" -max_memory %d -precision %d -verbose n -c %s -m %d -n %d -w R -a "2^%d" -b "2^%d" -d 1.0 -arg_type x -which R1_1,R1_2,R2_1,R2_2 -p %d %s data/pro_%s_R_small.txt', path, max_memory, precision, nice_number(c, 20), m, n, e, e, p, output_type, generate_name(c, m, n));
		end
		fprintf('%s\n', command);
		[ ...
		status, result ...
		] = system(command, '-echo');
		if (status == 1)
			return;
		end
	end
	files = {sprintf('pro_%s_S1.txt', generate_name(c, m, n)), sprintf('pro_%s_R.txt', generate_name(c, m, n)), sprintf('pro_%s_R_small.txt', generate_name(c, m, n))};
	zip(sprintf('data/pro_%s_%s.zip', generate_name(c, m, n), name), files, 'data');
	delete(sprintf('data/pro_%s_*.txt', generate_name(c, m, n)));
end
