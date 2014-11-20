#
# Copyright (c) 2014, Ross Adelman, Nail A. Gumerov, and Ramani Duraiswami
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#

#
# This script mass-launches several MATLAB instances, and in each of them, runs
# sandbox.m and run_tasks.m.  The calling convention for this script is:
#
# perl run_tasks.pl start_at num_processes skip_every
#
# For example,
#
# perl run_tasks.pl 1 8 8
#
# would launch eight MATLAB instances and run the following eight lines of
# code, one in each instance:
#
# run_tasks(1, 8, true);
# run_tasks(2, 8, true);
# run_tasks(3, 8, true);
# run_tasks(4, 8, true);
# run_tasks(5, 8, true);
# run_tasks(6, 8, true);
# run_tasks(7, 8, true);
# run_tasks(8, 8, true);
#
# Unlike when you run run_tasks by hand, the value for exit_after is always
# true, so MATLAB will exit after run_tasks is done.  Now, suppose you had
# eight machines, each with eight cores.  You could call run_tasks.pl on each
# of those machines:
#
# perl run_tasks.pl 1 8 64
# perl run_tasks.pl 9 8 64
# perl run_tasks.pl 17 8 64
# perl run_tasks.pl 25 8 64
# perl run_tasks.pl 33 8 64
# perl run_tasks.pl 41 8 64
# perl run_tasks.pl 49 8 64
# perl run_tasks.pl 57 8 64
#
# In this example, there'd be 64 MATLAB instances running.  Thus, this script
# provides some simple, course-grained parallelization.
#
$argc = $#ARGV + 1;
if ($argc != 3)
{
	print("bad calling convention...\n");
	exit(0);
}
$start_at = eval($ARGV[0]);
$num_processes = $ARGV[1];
$skip_every = $ARGV[2];
for $i (1 .. $num_processes)
{
	system("matlab -nosplash -nodesktop -nodisplay -r \"sandbox(); run_tasks(" . ($start_at + $i - 1) . ", " . $skip_every . ", true);\" > /dev/null &");
}
system("top");
