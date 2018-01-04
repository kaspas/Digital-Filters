clear;
close;

%% Initialize
%white noise
sigmav=0.54;
n=2^15  ;
v=sqrt(sigmav)*rand(n,1);
v=v-mean(v);

%create periodical interference

fo=1/4;
phi=pi/2;
A=4.2;
i=1:n;
x(i,1)=A*(sin(2*pi*fo*i+phi)+cos(4*pi*fo*i+phi)+cos(7*pi*i+phi/3));
s=x+v;

delta=8;%delay
M=100;%number of filter coefficients 

%u(n)=s(n-delta)
u=zeros(size(s));
u(delta+1:end)=s(1:end-delta);
%% Calculate correlations
[r,lags]=xcorr(u,u,M+1,'unbiased');
r=r(lags>=0);
R=toeplitz(r(1:end-1));
rbar=r(2:end);%Autocorrelation vecPor without r(0)

%%Wiener filter Coefficients

wo=R\rbar;
%% Call levinson function and compare with matlab functions
[a G L P]=my_levinson(r,M);
[matlab_a, matlab_P,matlab_G]=levinson(r,M);
fprintf('Filter coefficient MSE: %e \n', immse(a', matlab_a));
fprintf('Reflection coefficient MSE:  %e \n', immse(G, matlab_G));
fprintf('Prediction Error power MSE: %e \n', (P(end) - matlab_P) ^ 2);

%% Orthogonalization of u(n).
% Gram-Schmidt orthogonalization algorithm. 
b = zeros(n, M +1);
b(1:M) = u(1:M);

for i = M +1:n 
    b(i,:) = L * u(i:- 1:i - M);
end

% Autocorrelation matrix of backward errors.
D=diag(P);

p = zeros(M + 1, 1);

% Cross-correlation between each signal b[i] and desired signal s
for i = 1:M + 1
    [rb, lags] = xcorr(b(:, i), s, M, 'unbiased');
    p(i) = rb(lags == 0);
end

go = D\p; % The optimal parameters for the joint process estimator.
fprintf('L^T*go = wo MSE: %e \n', immse(L.' * go, wo));

%% Filters.
% Wiener.
figure(1);

y = zeros(n, 1);
for i = M + 2:n
    y(i) = wo' * u(i - 1:- 1:i - M - 1);
end
y = y(1:length(x));
e = s - y;
subplot(2,1,1);
plot((e-v).^2);
title('Error estimation Wiener');
xlabel('steps');
ylabel('(e-v)^2');
% JPE.
y = zeros(n, 1);
y(1:M) = u(1:M);
for i = M + 1:n
    y(i) = b(i,:) *go;
end
subplot(2,1,2);
e = s - y;
plot((e-v).^2);
title('Error estimation JPE');
xlabel('steps')
ylabel('(e-v)^2');

    
