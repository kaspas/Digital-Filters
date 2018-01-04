clear all;
s=matfile('sound.mat');
noise=matfile('noise.mat');
d=s.d;
n=length(d);
fs=s.fs;
u=noise.u;
sigmav=0.72;
%% filter
nCoeff=3;
a = xcorr(u,u,nCoeff-1,'unbiased');
a = a(nCoeff:(2*nCoeff-1));
R = toeplitz(a);
p=[sigmav;0;0];
% [D V]=eig(R);
% V=V*[1;1;1];
% 2/max(V)
wo=mldivide(R,p);
y=zeros(n,1);  
for i=3:n
  y(i) = u(i:-1:i-2)' * wo; % filter
end
%% clear sound
e=zeros(n,1);
e= d-y;
sound(e,fs);
