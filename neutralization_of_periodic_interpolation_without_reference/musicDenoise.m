clear; close;
M=130;
delta=100;
load music.mat
n=length(s);
d=s;
u=zeros(size(d));
u(delta+1:end)=s(1:end-delta);
[r,lags]=xcorr(u,u,M+1,'unbiased');
r=r(lags>=0);

[a G L P]=my_levinson(r,M);

%% Cross correlation.
% Calculate the cross correlation between each signal b[i] and the desired
% signal s. Because b's size, n x (M+1), is too large (>4GB in current
% example) we caclulate it seperately for each coefficient and don't save
% the values afterwards.[a, G, L, P] = my_levinson(r,M);
p=zeros(M+1,1);
for i=1:M+1
    b=filter(L(i,1:i),1,u);
    b(1:i)=u(1:i);
    [rb,lags]=xcorr(b,s,1,'unbiased');
    p(i)=rb(lags==0);
end
D=diag(P);
go=D\p;
%% Clear sound
y=zeros(n,1);
for i=1:M+1
    b=filter(L(i,1:i),1,u);
    b(1:i)=u(1:i);
    y=y+go(i)*b;
end
e = s - y;
%% start sound
sound(e, fs);
