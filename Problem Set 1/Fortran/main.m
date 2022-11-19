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
    !./opt_vfi
elseif strcmp(os,'win')==1
    !opt_vfi.exe
else
    error('Specify which operating system you are using; variable "os"')
end

%End timer
toc

load k_grid.txt
load valuefunction.txt
load k_policy.txt

figure
plot(k_grid, valuefunction)
xlabel('k')
ylabel('V(k)')
title('Value Function')

figure
plot(k_grid, k_policy)
xlabel('k')
ylabel('g(k)')
title('Policy Function')

