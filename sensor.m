function [f,g,H,J] = sensor(x,params)
% x: value of the optimization variable
% params: structure containing the parameters of the objective function
% params.e = light intensity values
% params.r = sensor resistance values
% f,g,H,J: objective function, its gradient, hessian and Jacobian
a = x(1);
b = x(2);
R = params.r;
E = params.e;

f = 0.5*sum((R - b.*E.^(-a)).^2);

g1 = sum((R - b.*E.^(-a)).*(b.*log(E).*E.^(-a)));
g2 = sum((R - b.*E.^(-a)).*(-1).*E.^(-a));
g = [g1;g2];

H11 = sum(b*R.*(log(E)).^2.*E.^(-a));
H12 = sum(-b*log(E).*E.^(-2*a)+(R-b./E.^a).*(log(E)./E.^a)  );
H21 = H12;
H22 = sum( E.^(-2*a) );
H = [H11 H12;H21 H22];

J1 = b*log(E)./E.^a;
J2 = -E.^(-a);
J = [J1', J2'];
end