%% Init
clear all;

load speakerA;
fsA=fs;
load speakerB;

fs=fsA;
clear fsA;

M=6600;
tic
[ylms,wlms,elms]=my_lms(u,d,M);
toc
tic
[ynlms,wnlms,enlms]=my_nlms(u,d,M);
toc


%[yrls,wrls,erls]=my_rls(u,d,M);

% fprintf('Playing sound resulting from RLS.\nPress Enter to stop sound.\n');
 
% sound(erls,fs);
 
 fprintf('Playing sound resulting from NLMS.\nPress Enter to stop sound.\n');
 
 sound(enlms,fs);