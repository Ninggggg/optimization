function [f,g,H,J] = rosenbrock(x, params)
    n = params.n;
    b = params.b;
    for i=1:n-1
        f = b*(x(i+1)-x(i)^2)^2 + (1-x(i))^2;
    end
    
    if n==2
        g = [-4*b*x(1)*(x(2)-x(1)^2)-2*(1-x(1)); 2*b*(x(2)-x(1)^2)];
        H = [12*b*x(1)^2-4*b*x(2)+2, -4*b*x(1); -4*b*x(1), 2*b];
        J = [sqrt(2*b)*(-2)*x(1), sqrt(2*b); -sqrt(2), 0];
%     elseif n>=3
%         g = zeros(n,1);
%         H = zeros(n, n);
%         g(1) = -4*b*x(2)*( x(2)-x(1)^2)-2*(1-x(1));
%         % H(1,1) = 
%         for i=2:n-1
%             g(i) = -4*b*x(i)*(x(i)-x(i-1)^2)-2*(1-x(i-1))+2*b*(x(i+1)-x(i)^2);
%         end
%         g(n) = 2*b*(x(n)-x(n-1)^2);
     end
end
    