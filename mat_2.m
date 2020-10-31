clear all
close all

E = [5013 2415 1558 1000 820 621 433 201 105 55];
R = 0.001*[0.141 0.329 0.525 0.970 1.140 1.511 2.362 5.224 12.826 25.512];
params.e = E;
params.r = R;
x0 = [1, 2]';
x1 = [3, 3]';
options.method = 'gradient';
options.step = 'variable';
options.const = 1e-3;
options.beta = 0.75;
options.TolX = 1e-8;
options.TolF = 1e-8;
options.TolG = 1e-8;
options.MaxIter = 1e4;

[xh, result, xval] = optimdescent('sensor',params,options,x0);
options.method = 'conjugate';
[xh2, result2, xval2] = optimdescent('sensor',params,options,x0);
options.method = 'Newton';
[xh3, result3, xval3] = optimdescent('sensor',params,options,x0);
options.method = 'Quasi-Newton';
[xh4, result4, xval4] = optimdescent('sensor',params,options,x0);
options.method = 'Guass-Newton';
[xh5, result5, xval5] = optimdescent('sensor',params,options,x0);
options.method = 'Levenberg-Marqardt';
[xh6, result6, xval6] = optimdescent('sensor',params,options,x0);

figure
hold on;plot(log(result.crit));plot(log(result2.crit));plot(log(result3.crit));
plot(log(result4.crit));plot(log(result5.crit));plot(log(result6.crit));hold off;
legend('steepest gradient','conjugate','Newton','Quasi-Newton','Guass-Newton','Levenberg-Marqardt')
title('value of the criterion with all methods')

figure
hold on;plot(log(result.grad));plot(log(result2.grad));plot(log(result3.grad));
plot(log(result4.grad));plot(log(result5.grad));plot(log(result6.grad));hold off;
legend('steepest gradient','conjugate','Newton','Quasi-Newton','Guass-Newton','Levenberg-Marqardt')
title('norm of the gradient with all methods')

disp(['The computation time of steepest gradient is ',num2str(result.time)])
disp(['The computation time of conjugate gradient is ',num2str(result2.time)])
disp(['The computation time of Newton method is ',num2str(result3.time)])
disp(['The computation time of Quasi-Newton method is ',num2str(result4.time)])
disp(['The computation time of Guass-Newton method is ',num2str(result5.time)])
disp(['The computation time of Levenberg-Marqardt method is ',num2str(result6.time)])

x1 = [0:0.01:2]';
x2 = [0:0.01:4]';
[X1,X2] = meshgrid(x1,x2);
size1 = length(x1);
size2 = length(x2);
F = zeros(size1,size2);
for i = 1:size1
    for j =1:size2
        [f,g,H] = sensor([x1(i);x2(j)],params);
        F(i,j)=f;
    end
end
figure
subplot(2,3,1);contour(X1',X2',log10(F),500);hold on;
plot(xval(1,:),xval(2,:),'-g','LineWidth',2);title('gradient descent');
text(x0(1),x0(2),['point initial (', num2str(x0(1)),',',num2str(x0(2)),')'],'HorizontalAlignment','right');
text(xh(1),xh(2),['point final (', num2str(xh(1)),',',num2str(xh(2)),')'],'HorizontalAlignment','left');hold off;


subplot(2,3,2);hold on;contour(X1',X2',log10(F),500);
plot(xval2(1,:),xval2(2,:),'-g','LineWidth',2);title('conjugate descent');
text(x0(1),x0(2),['point initial (', num2str(x0(1)),',',num2str(x0(2)),')'],'HorizontalAlignment','right');
text(xh2(1),xh2(2),['point final (', num2str(xh2(1)),',',num2str(xh2(2)),')'],'HorizontalAlignment','right');hold off;

subplot(2,3,3);hold on;contour(X1',X2',log10(F),500);
plot(xval3(1,:),xval3(2,:),'-g','LineWidth',2);title('Newton method');
text(x0(1),x0(2),['point initial (', num2str(x0(1)),',',num2str(x0(2)),')'],'HorizontalAlignment','right');
text(xh3(1),xh3(2),['point final (', num2str(xh3(1)),',',num2str(xh3(2)),')'],'HorizontalAlignment','right');hold off;

subplot(2,3,4);hold on;contour(X1',X2',log10(F),500);
plot(xval4(1,:),xval4(2,:),'-g','LineWidth',2);title('Quasi-Newton');
text(x0(1),x0(2),['point initial (', num2str(x0(1)),',',num2str(x0(2)),')'],'HorizontalAlignment','right');
text(xh4(1),xh4(2),['point final (', num2str(xh4(1)),',',num2str(xh4(2)),')'],'HorizontalAlignment','right');hold off;

subplot(2,3,5);hold on;contour(X1',X2',log10(F),500);
plot(xval5(1,:),xval5(2,:),'-g','LineWidth',2);title('Guass-Newton');
text(x0(1),x0(2),['point initial (', num2str(x0(1)),',',num2str(x0(2)),')'],'HorizontalAlignment','right');
text(xh5(1),xh5(2),['point final (', num2str(xh5(1)),',',num2str(xh5(2)),')'],'HorizontalAlignment','right');hold off;

subplot(2,3,6);hold on;contour(X1',X2',log10(F),500);
plot(xval6(1,:),xval6(2,:),'-g','LineWidth',2);title('Levenberg-Marqardt');
text(x0(1),x0(2),['point initial (', num2str(x0(1)),',',num2str(x0(2)),')'],'HorizontalAlignment','right');
text(xh6(1),xh6(2),['point final (', num2str(xh6(1)),',',num2str(xh6(2)),')'],'HorizontalAlignment','right');hold off;