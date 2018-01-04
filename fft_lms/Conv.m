clear all
close all

n1 = 1024; 
n2 = 1024;
% input 
x = randn(n1,1) + 1i*randn(n1,1);
y = randn(n2,1) + 1i*randn(n2,1);

%% Convolution

Co1=conv(x,y);

%% Toeplitz convolution
Y=toeplitz([y;zeros(n2-1,1)],[y(1);zeros(n2-1,1)]);
Co2=Y*x;
fprintf('Toeplitz convolution error : %e\n', norm(Co1 - Co2))

%% Circulant convolution
yr=[y(1) zeros(1,n2-1) y(end:-1:2).'];
C=toeplitz([yr(1) fliplr(yr(2:end))], yr );
Co3=C*[x ; zeros(n1-1,1)];
fprintf('Circulant convolution error : %e\n', norm(Co1 - Co3))

%% Inverse FFT/FFT convolution
yp=[y; zeros(n2-1,1)];
xp=[x; zeros(n1-1,1)];
Co4=ifft(fft(yp).*fft(xp));
fprintf('Inverse FFT/FFT convolution error : %e\n', norm(Co1 - Co4))