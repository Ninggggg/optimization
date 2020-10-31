function [xh,result,xval] = optimdescent(critfun,params,options,x0)
% critfun : name of fonction .m evaluating
% the objective function, its gradient and hessien
% params : parameter for objective function evaluation
% options : necessary options for the algorithms
% options.method : direction type, ¡¯gradient¡¯, ¡¯Newton¡¯, ¡¯Quasi-Newton¡¯, ...
% options.step : ¡¯fixed¡¯, ¡¯variable¡¯
% options.const : Armijo constant
% options.beta : backtracking rate
% options.TolX, options.TolF, options.TolG, options.MaxIter : stopping thresholds
% x0 : initial point
% xh : final point
% result : structure containing the results
% result.iter : number of iterations performed
% result.crit : value of the criterion by iteration
% result.grad : norm of the gradient by iteration
% result.temp : calculation time
% result.stop : stop condition(¡¯TolX¡¯, ¡¯TolG¡¯, ¡¯TolF¡¯, ¡¯Maxiter¡¯)
% xval : values ??of iterates
tic
iterate = 1;

eta_k = 0.01;

x_pre = x0;
[fk, gk, Hk, Jk] = feval(critfun, x0, params);
f_pre = fk;     % f0
g_pre = gk;     % g0
dk = -g_pre;    
alpha = params.alpha;
xk = x0;
% fk_ = feval(critfun, xk+alpha*dk, params);
% while fk_ >= (fk + options.const*alpha*gk'*dk)
%     alpha = options.beta*alpha;
%     fk_ = feval(critfun, xk+alpha*dk, params);
% end
result.iter = 1;

xval(:,result.iter) = x0;
xk = x0 + alpha*dk;
Bk_inv = eye(length(x0));
result.crit(:,result.iter) = fk;
result.grad(:,result.iter) = norm(gk);

while iterate
    
[fk, gk, Hk, Jk] = feval(critfun, xk, params);
dk = -gk;

switch options.method
    case 'gradient'
        dk = -gk;
    case 'conjugate'  % Fletcher-Reeves version 
        options.beta = norm(gk)/norm(g_pre);
        dk = -gk + options.beta*dk;
    case 'Newton'
        dk = -inv(Hk)*gk;
    case 'Quasi-Newton'
        sk = xk - x_pre;
        yk = gk - g_pre;
        rho_k = 1/(yk.'*sk);
        Vk = eye(length(x0)) - rho_k*yk*sk.';
        Bk_inv = Vk.'*Bk_inv*Vk + rho_k*sk*sk.';
        dk = -Bk_inv*gk;
    case 'Guass-Newton'
        dk = -inv(Jk'*Jk)*gk;
    case 'Levenberg-Marqardt'
        lambda = eta_k*norm(gk);
        dk = -inv(Jk'*Jk + lambda*diag(diag(Jk'*Jk)))*gk;
end

if dk'*gk > 0
    dk = -dk;
end

% options.beta = norm(gk)/norm(g_pre);
switch options.step
    case 'fixed'
        alpha = params.alpha;
    case 'variable'
        alpha = 1;
        % fk_ = rosenbrock(xk+alpha*dk, params); % calculate the f value every loop
        fk_ = feval(critfun, xk+alpha*dk, params);
        while fk_ >= fk + options.const*alpha*gk'*dk
            alpha = options.beta*alpha;
            fk_ = feval(critfun, xk+alpha*dk, params);
        end
end

x_pre = xk;
xk = xk + alpha * dk;
result.iter = result.iter + 1;

if  norm(gk) < options.TolG
    iterate = 0;
    result.stop = 'TolG';
end

if  result.iter > options.MaxIter
    iterate = 0;
    result.stop = 'MaxIter';
end

if abs((fk-f_pre)/f_pre) < options.TolF
    iterate = 0;
    result.stop = 'TolF';
end
    
if norm(xk-x_pre)/norm(x_pre) < options.TolX
    iterate = 0;
    result.stop = 'TolX';
end

result.crit(:,result.iter) = fk;
result.grad(:,result.iter) = norm(gk);
xval(:,result.iter) = xk;
f_pre = fk;
g_pre = gk;

end
xh = xk;
result.time = toc;