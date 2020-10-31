function [f,g,H,J] = fgh_log(S,params)
beta = params.beta;
J = -sqrt(2)*params.K;

if isempty(find(S<0))
   f = (params.y-params.K*S)'*(params.y-params.K*S) + params.beta*(params.D*S)'*(params.D*S) - params.mu*sum(log(S));
else
    f = +Inf;
end
g = -2 * params.K' * (params.y-params.K*S) + 2*beta*params.D'*(params.D*S) - params.mu./S;
H = 2*params.K'*params.K + 2*beta*params.D'*params.D + diag(params.mu./(S.^2));
end