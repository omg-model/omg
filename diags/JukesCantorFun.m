function [F,jacF] = JukesCantorFun(t,p)

% F = f(t)
F = 0.5.*(1-exp(t.*p));

if nargout > 1 % use Jacobian
    % first derivative of F wrt t
    jacF = -0.5.*p.*exp(t.*p);
end


