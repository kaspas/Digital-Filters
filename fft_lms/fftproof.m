% "Proof" by MATLAB
% A simple technique to develop and verify the steps of a proof 
% using random data input
%
% N P P
% Cornell U 
% Sept 1992
%

clear

n = 8; % any even

% input 
x = randn(n,1) + 1i*randn(n,1);
% correct answer

ys = fft(x);

% root of unity
w = @(n,e) exp(-2*pi*1i.*e/n);

k = (0:n-1)';

% DFT proof steps
y = zeros(n,1);

for j = 0:n-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end

fprintf('DFT : %e\n', norm(y - ys))

% split output top bottom
y = zeros(n,1);
for j = 0:n/2-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end
for j = n/2:n-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end
fprintf('split output top bottom : %e\n', norm(y - ys))

% split input even odd
y = zeros(n,1);

%% fft Proof
fy1=zeros(n/2-1,1);
fy2=zeros(n/2-1,1);
k = (0:n/2-1)';
for j=0:n/2-1
  fy1(j+1)=sum(w(n/2,k*(j)).*x(2*k+1));
  fy2(j+1)=sum(w(n/2,k*(j)).*x(2*k+2));
end
W=diag(w(n,k));
fy=[fy1+W*fy2;fy1-W*fy2];
fprintf('FFT algorithm : %e\n',norm(fy-ys));

%%
k = (0:n/2-1)';

for j = 0:n/2-1
  y(j +1) = sum(w(n,j*2*k) .* x(2*k +1)) + sum(w(n,j*(2*k+1)) .* x(2*k+1 +1));
end
for j = n/2:n-1
  y(j +1) = sum(w(n,j*2*k) .* x(2*k +1)) + sum(w(n,j*(2*k+1)) .* x(2*k+1 +1));
end

fprintf('split input even odd : %e\n', norm(y - ys))
%% Recursive fft
tic;
ry=r_fft(x);
toc;
fprintf('recursive fft : %e\n',norm(ry-ys));
%%  
% apply w identities
% etc
% ...
% ...

% to complete the proof
fe = fft(x((0:2:n-1) +1));
fo = fft(x((1:2:n-1) +1));  
wfo = w(n,(0:n/2-1)') .* fo; 
y = [fe + wfo; fe - wfo];
fprintf('done my proof : %e\n',norm(fy-y));
fprintf('done : %e\n', norm(y - ys))

