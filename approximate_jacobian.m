%Implementation of finite difference approximation
%for Jacobian of multidimensional function
%INPUTS:
%fun: the mathetmatical function we want to differentiate
%X: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%J: approximation of Jacobian of fun at x
function [J,num_eval] = approximate_jacobian(fun,X)
    m = length(fun(X));
    n = length(X);
    J = zeros(m,n);
    %set the step size to be tiny
    delta = 1e-6;
    num_eval = 1;
    for col=1:n
        mask = zeros(n,1);
        mask(col) = delta;
        %compute the function at different points near x
        f_left = fun(X-mask);
        f_right = fun(X+mask);
        %approximate the first derivative
        J(:,col) = (f_right-f_left)/(2*delta);
        num_eval = num_eval + 2;
    end
    
end