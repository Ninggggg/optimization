close all
clear all
clc

E = [5013 2415 1558 1000 820 621 433 201 105 55];
R = 0.001*[0.141 0.329 0.525 0.970 1.140 1.511 2.362 5.224 12.826 25.512];
N = length(E);

% x = [-2, 1]';
params.e = E;
params.r = R;

a = [0:0.01:2]';
b = [0:0.01:4]';

N1 = length(a);
N2 = length(b);
F = zeros(N1, N2);

for i=1:N1
    for j=1:N2
        x = [a(i), b(j)]';
        F(i,j) = sensor(x, params);
    end
end

[X1,X2] = meshgrid(a,b);
contour(X1',X2',log10(F),500)
figure
meshc(X1',X2',log10(F))