function X = r_fft(x)
  n=length(x);
  w = @(n,e) exp(-2*pi*1i.*e/n);    

%only works if N = 2^k
%x1=even , x2=odd
  x1=x(1:2:end);  %even
  x2=x(2:2:end);  %odd
  if n==1 
    X=x;
  else
    X1=r_fft(x1);
    X2=r_fft(x2);
    X=zeros(n,1);
    W=w(n,(0:(n/2)-1)');
    z=W .* X2;
    X=[(X1+z);(X1-z)];
  end
end