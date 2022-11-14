%This m-file codes the Problem Set 1 using a Fortran .exe

%Preliminaries:
clear
close all
clc

%Start timer
tic

%Determine which operating system you're using:
findos %returns os='win' or os='unix'


%Call the executable:
if strcmp(os,'unix') == 1
    !./vfi
elseif strcmp(os,'win')==1
    !vfi.exe
else
    error('Specify which operating system you are using; variable "os"')
end

%End timer
toc

