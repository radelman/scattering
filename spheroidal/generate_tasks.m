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
% This will generate 2 * M * N tasks in tasks/ready.  The precisions and
% numbers of coefficients were chosen for c = 10 and c = 25.  These will also
% work nicely for any values of c below that.  Above that, you may need to
% increase them.
%
% For example, suppose you want to compute the spheroidal wave functions for
% c = 15.0, m = 0, 1, ..., 4, and n = m, m + 1, ... m + 9.  You can do:
%
% generate_tasks(15.0, 5, 10, 0);
%
% This will generate 100 tasks in tasks/ready, 50 for the prolate case and 50
% for the oblate case.
%
% As another example, supose you want to do the same thing, but for both c =
% 15.0 and c = 17.0.  You can do:
%
% generate_tasks(15.0, 5, 10, 0);
% generate_tasks(17.0, 5, 10, 100);
%
% Note how, in the second line, the last argument is 100.  The tasks created
% in this call will be numbered 101, 102, and so on.  This is good because the
% last task created by the preceding line was 100.  In other words, task
% numbers will start at task_number_start + 1.  This allows you to generate
% multiple sets of tasks without overwriting previously generated ones.
%
function generate_tasks(c, M, N, task_number_start)
	task_number = task_number_start;
	for m = 0 : M - 1
		for n = m : m + N - 1
			task_number = task_number + 1;
			fid = fopen(sprintf('tasks/ready/task%05d.m', task_number), 'w');
			fprintf(fid, 'pro_calculate_everything(sprintf(''%%s/sphwv'', pwd()), 2000, 500, %s, %d, %d, 10, -300, 10, -300, 10, -300);\n', nice_number(c, 20), m, n);
			fprintf(fid, 'pro_calculate_functions('''', sprintf(''%%s/sphwv'', pwd()), 2000, 300, %s, %d, %d, 1.0 / 2048.0, ''theta/pi'', -9, 9.0, 1.0 / 256.0, ''x'', 17);\n', nice_number(c, 20), m, n);
			fprintf(fid, 'everything = pro_load_everything(''data'', %s, %d, %d, '''');\n', nice_number(c, 20), m, n);
			fprintf(fid, 'pro_save_something(''saved'', everything);\n');
			fclose(fid);
			task_number = task_number + 1;
			fid = fopen(sprintf('tasks/ready/task%05d.m', task_number), 'w');
			fprintf(fid, 'obl_calculate_everything(sprintf(''%%s/sphwv'', pwd()), 2000, 500, %s, %d, %d, 10, -300, 10, -300, 10, -300, 10, -300);\n', nice_number(c, 20), m, n);
			fprintf(fid, 'obl_calculate_functions('''', sprintf(''%%s/sphwv'', pwd()), 2000, 300, %s, %d, %d, 1.0 / 2048.0, ''theta/pi'', 8.0, 1.0 / 256.0, ''z'', 17);\n', nice_number(c, 20), m, n);
			fprintf(fid, 'everything = obl_load_everything(''data'', %s, %d, %d, '''');\n', nice_number(c, 20), m, n);
			fprintf(fid, 'obl_save_something(''saved'', everything);\n');
			fclose(fid);
		end
	end
	fprintf('There are now %d tasks...\n', task_number);
end
