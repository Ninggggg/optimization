clear all
close all
%
N  = 200;
M  = 5000;
tmax = 10;
Tmax = tmax/4;

t    =  linspace(0,tmax,M);
T    =  linspace(0,Tmax,N+1);
T(1) = []; %remove first time instant to avoid divition over 0

stdT = Tmax/30;
muT = Tmax/4;

rho = 0.02;
S = exp(-abs(((T-muT)/stdT).^2));

S = S(:);

[K] = simrelax(T,t,'T2');

y = K*S;

SNR = 100;
y = y + (max(y)/SNR)*randn(size(y));

figure(1);
subplot(121)
plot(T,S)
xlabel('Relaxation time')

subplot(122),
plot(t,y)

%% log barrier method
tic
options.method = 'Newton';
options.step = 'variable';
options.const = 1e-3;
options.beta = 0.75;
options.TolX = 1e-8;
options.TolF = 1e-8;
options.TolG = 1e-8;
options.MaxIter = 1e4;

e = ones(N-1,1);
D = spdiags([-e e],0:1,N-1,N);
D = full(D);

mu0 = 1;

s0 = 1/N*ones(N,1);
s_est_pre = s0;

params.D = D;
params.mu = mu0;
params.K = K;
params.y = y;
params.beta = 1e-3;

flag = 1;
tol = 1e-2;
s_est = optimdescent('fgh_log',params,options,s0);    
figure;plot(T,s_est);
count = 0;
figure
hold on
while flag
    s_est = optimdescent('fgh_log',params,options,s0);    
    if norm(s_est-s_est_pre) < tol
        flag = 0; 
    end
    plot(T,s_est)
    s_est_pre = s_est;
    params.mu = params.mu/10;
    count = count + 1;
end
toc
hold off
figure
subplot 211;plot(T,S);title('Real S')
subplot 212;plot(T,s_est);title('estimated by log barrier method')
