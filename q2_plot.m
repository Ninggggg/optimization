clear all
close all

E = [5013 2415 1558 1000 820 621 433 201 105 55];
R = 0.001*[0.141 0.329 0.525 0.970 1.140 1.511 2.362 5.224 12.826 25.512];

C = [E',ones(10,1)];
d = R;
[xh,xval,residual,exitflag] = lsqlin(C,d);

plot(E,R,'*');
hold on;
alpha = 1.1477;
beta = 2.5508;

E1 = 50:1:6000+50;
R1 = beta*E1.^(-alpha);
plot(E1,R1)