clear all
close all
clc

x0 = [-3, 10]';
x1 = [2, 1]';
params.n = 2;
params.b = 2;
options.method = 'gradient';
options.step = 'variable';
options.const = 1;
options.beta = 0.75;
options.TolX = 1e-8;
options.TolF = 1e-8;
options.TolG = 1e-8;
options.MaxIter = 1e4;
options.const = 1e-3;


[xh, result, xval] = optimdescent('rosenbrock',params,options,x0);

options.method = 'conjugate';
[xh2, result2, xval2] = optimdescent('rosenbrock',params,options,x0);

options.method = 'Newton';
[xh3, result3, xval3] = optimdescent('rosenbrock',params,options,x0);

options.method = 'Quasi-Newton';
[xh4, result4, xval4] = optimdescent('rosenbrock',params,options,x0);

options.method = 'Guass-Newton';
[xh5, result5, xval5] = optimdescent('rosenbrock',params,options,x0);

options.method = 'Levenberg-Marqardt';
[xh6, result6, xval6] = optimdescent('rosenbrock',params,options,x0);

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