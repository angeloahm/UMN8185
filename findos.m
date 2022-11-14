%This file finds the OS of the computer on which the program runs.
%This is needed for the usage of the FORTRAN executables.

os_raw = computer('arch');

switch os_raw
     case 'win32'
		os = 'win';
     case 'win64'
		os = 'win';
     case 'glnx86'
		os = 'unix';
     case 'glnxa64'
		os = 'unix';
     case 'maci64'
		os = 'unix';
end
