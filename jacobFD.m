function [J] = jacobFD(f,x,delx)
% Calculates the Jacobian of the
% system of non-linear equations:
% f(x) = 0, through finite differences.
% The Jacobian is built by columns
[m, ~] = size(x);
J = zeros(m,m);
for j = 1:m
    xp = x;
    xm = x;
    xp(j) = x(j) + delx;
    xm(j) = x(j) - delx;
    J(:,j) = (feval(f,xp)-feval(f,xm))/delx/2;
end; 
% end function