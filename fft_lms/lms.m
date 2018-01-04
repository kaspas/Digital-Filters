clear all
close all

n = 2^15;
trials = 1;
sigma = 0.57;

v  = zeros(trials,  n);
for t = 1:trials
    v(t,:) = sqrt(sigma)*randn(n,1); v(t,:) = v(t,:) - mean(v(t,:));
end

uin = zeros(trials, n);
a = 0.34;
for t = 1:trials
   
    uin(t,1) = v(t, 1);
    for i=2:n
        uin(t,i) = -a * uin(t, i-1) + v(t, i);
    end
    
end


% generate the system
d = plant(uin);


%% adaptation
mu =0.6*0.0005;
M = 2^10; % taps


J = zeros(n, 1);
Jm = J; Js = J; Js2=J;
for t=1:trials
  
  tic  
  %% loop version
  y = zeros(n,1);
  e = zeros(n,1);
  u = uin(t,:)'; % input to the system
  w = zeros(M, 1);
  
  for i = 2*M:M:n
    g = zeros(M,1);
    for j=M-1:-1:0
      r = i-j;
      
      y(r) = w'*u((r-M+1):r);
      e(r) = d(t, r) - y(r);
      
      g = g + u((r-M+1):(r))*e(r);
    end
    w = w + mu*g;
    
  end
  J = J + e.^2;
  
  
  toc
tic
  %% matrix version
  ym = zeros(n,1);
  em = zeros(n,1);
  wm = zeros(M,1);
  for i = 2*M:M:n
    
    S = toeplitz(u((i-M+1):i), u((i-M+1):-1:(i-2*M+2)));
    ym(i-M+1:i) = S * wm;
    em(i-M+1:i) = d(t, i-M+1:i) - ym(i-M+1:i)';
    
    g = S' * em(i-M+1:i);
    wm = wm + mu * g;
  end
  Jm = Jm + em.^2;
  
  toc
  tic
  %% FFT version 1 contstrained
  ys = zeros(n,1);
  es = zeros(n,1);
  Ws = zeros(2*M,1);

 for i = 2*M:M:n
    
    U = fft(u(i-2*M+1:i));
    C = ifft(Ws .* U);
    ys(i-M+1:i) = C(M+1:end);
    es(i-M+1:i) = d(t, i-M+1:i) - ys(i-M+1:i)';
    g = ifft(fft([zeros(M,1);es(i-M+1:i)]) .* conj(U));
    G=fft([g(1:end-M);zeros(M,1)]);
    Ws = Ws + mu * G;

  end
  
  Js = Js + es.^2;
toc
tic
  %% FFT version 2 unconstrained(same mu)
  ys = zeros(n,1);
  es = zeros(n,1);
  Ws = zeros(2*M,1);
 for i = 2*M:M:n
    U=(fft(u(i-2*M+1:i)));
    C=ifft(Ws.* U);
    ys(i-M+1:i)=C(M+1:end);
    es(i-M+1:i) = d(t, i-M+1:i) - ys(i-M+1:i)';
    g=fft([zeros(M,1);es(i-M+1:i)]).*conj(U);
    Ws=Ws+mu*g;
  end
  
   Js2 = Js2 + es.^2;
   toc

end


J = J / trials;
Jm = Jm / trials;
Js = Js / trials;
Js2 = Js2 / trials;


% 
% 
% figure
% plot([J - Jm])
% xlabel('Time steps')
% ylabel('error');
% title('Block LMS and FAST LMS error difference')

figure(1)
subplot(2,2,1)
hold on
plot([J])
xlabel('Time steps')
ylabel('error');
title('Learning Curve(Nested for Loop)')
subplot(2,2,2)
plot([Jm])
xlabel('Time steps')
ylabel('error');
title('Learning Curve(Matrix Version)')
subplot(2,2,3)
plot([Js])
xlabel('Time steps')
ylabel('error');
title('Learning Curve(Constrained)')
subplot(2,2,4)
plot([Js2])
xlabel('Time steps')
ylabel('error');
title('Learning Curve(Unconstrained)')
figure(2)
plot([J Jm Js Js2])
xlabel('Time steps')
ylabel('error');
title('Block LMS and FAST LMS error difference')
legend('nestedLoop','Loop && Matrix','Fast LMS constrained','Fast LMS Unconstrained(mu=0.6*mu');
fprintf('J-Jm error : %e\n', norm(J - Jm))
fprintf('J-JS error : %e\n', norm(J - Js))
fprintf('J-JS2 error : %e\n', norm(J - Js2))