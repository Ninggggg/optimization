close all 
clear all

x = [-2, 1]';
params.n = 2;
params.b = 2;

% [f, g, H] = rosenbrock(x, params);

step = 0.01;
x1 = [-4:step:4]';
x2 = [-5:step:20]';
[X1, X2] = meshgrid(x1,x2);
N1 = length(x1);
N2 = length(x2);

F = zeros(N1, N2);
for i=1:N1
    for j=1:N2
        x = [x1(i), x2(j)]';
        F(i,j) = rosenbrock(x, params);
    end
end

[X1,X2] = meshgrid(x1,x2);
meshc(X1', X2', log10(F));
title('Image 3D of Rosenbrock in log')

L = 100;
figure
contour(X1', X2', F, L);
hold on
plot(X1(1),X2(1),'ro')
% plot(x1(1),X2(1),F(X1(1),F(X2(1)))'ro')