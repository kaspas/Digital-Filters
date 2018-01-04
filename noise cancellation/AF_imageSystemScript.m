clear all;
n=1024;
sv=0.32;   
for mu=1.0 : 1.5 :5.5
        figure(2*mu)
    %% signal
    A=sqrt(0.15)*randn(n,1);
    A=A-mean(A);
    fi=pi/6;
    x=zeros(n,1); % initialize
    N=1:n;
    for i=1:n
      x(i)=A(i)*sin((pi/8)*i+fi);
    end
    u=zeros(n,1);
    % noise
    v=sqrt(0.32)*randn(n,1);
    v=v-mean(v);
    % d(n) signal
    d=zeros(n,1);
    for i=1:n
    d(i)=x(i)+v(i);%add noise
    end
    % sensor value
    u(1)=v(1);
    u(2)=0.25*u(1)+v(2);
    for i=3:n
        u(i)=0.25*u(i-1)-0.12*u(i-2)+v(i);
    end
    subplot (3,1,1);
    plot([d x u])
    legend({'d(n)', 'x(n)', 'u(n)'})
    %% fir filter
    A=[1 -0.25 0.12;-0.25 1.12 0;0.12 -0.25 1];
    B=[sv 0 0];
    r=mldivide(A,B');
    p=[sv;0;0];
    R_u=[r(1) r(2) r(3);r(2) r(1) r(2);r(3) r(2) r(1)];
    wo=mldivide(R_u,p);
    Jmin=0.395 - p'*inv(R_u)*p;
    [D V]=eig(R_u);
    V=V*[1;1;1];
    2/max(V);
    w=[-1.5; -1.5 ; -1.5];
   % mu=2;
    wt=zeros([3,n]);
    wt(:,1)=w;
    y=zeros(n,1);
    for i=3:n
      w = w + mu*(p-R_u*w); % Adaptation steps
      wt(:,i) = w;
      y(i) = u(i:-1:i-2)' * w; % filter
    end
    e=zeros(n,1);
    e=d-y;
    subplot(3,1,2);
    plot([x e])
    legend({'x(n)', 'e(n)'})
    %parameter error
    subplot(3,1,3); 
    % figure(3)
    we = (wt - wo*ones(1,n)).^2;
    e = sqrt(sum(we));
    semilogy(e);
    xlabel('time step n');
    ylabel('Parameter error');
    title('Parameter error');

end
