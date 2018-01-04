
%% Echo cancellation

%% Initialization

n=2^20;
sigmav1=0.42;
sigmav2=0.72;
v1=sqrt(sigmav1)*rand(n,1);
v2=sqrt(sigmav2)*rand(n,1);

u=zeros(1,n);
d=zeros(1,n);
x=zeros(1,n);
s=zeros(1,n);

for i=4:1:n
    u(i)= -0.87*u(i-1)-0.22*u(i-2)-0.032*u(i-3)+v1(i);
    x(i)= -0.57*x(i-1)-0.16*x(i-2)-0.08*x(i-3)+v2(i);
    s(i)= -0.13*u(i)+0.67*u(i-1)-0.18*u(i-2)+0.39*u(i-3);
    d(i)=s(i)+x(i);
end
M=4;

%% Wiener coefficients

r=[1 0.87 0.22 0.032;0.87 1.22 0.032 0;0.22 0.902 1 0;
    0.032 0.22 0.87 1]\[sigmav1;0;0;0];
R=toeplitz(r);
p=[-0.13 0.67 -0.18 0.39;0.67 -0.31 0.39 0;
    -0.18 1.06 -0.13 0;0.39 -0.18 0.67 -0.13]*r;
wo=R\p;

[ylms,wlms,elms]=my_lms(u,d,M);
[ynlms,wnlms,enlms]=my_nlms(u,d,M);
[yrls,wrls,erls]=my_rls(u,d,M);

Jlms = elms .^ 2;
Jnlms = enlms .^2;
Jrls= erls .^2;
figure(1);
subplot(3,1,1);
plot(Jlms);
xlabel('time steps');
ylabel('e.^2');
title('LMS')
subplot(3,1,2);
plot(Jnlms);
xlabel('time steps');
ylabel('e.^2');
title('NLMS');
subplot(3,1,3);
plot(Jrls);
xlabel('time steps');
ylabel('e.^2');
title('RLS');
fprintf('MSE wo,wlms: %e \n',immse(wo,wlms))
fprintf('MSE wo,wnlms: %e \n',immse(wo,wnlms))
fprintf('MSE wo,wrls: %e \n',immse(wo,wrls))

fprintf('MSE wlms,wnlms: %e \n',immse(wlms,wnlms))
fprintf('MSE wlms,wrls: %e \n',immse(wlms,wrls))
fprintf('MSE wnlms,wrls: %e \n',immse(wnlms,wrls))

